"""
周期性 NMSE/SNR 起伏源于训练阶段仅在滑窗片段上重置隐状态，而推理时整段前向导致隐状态漂移
新增 reconstruct_sequence_with_overlap 通过 Hann 加权的滑窗推理复现训练分布以压制周期误差。
评估循环复用训练集的窗口长度与步长，逐段重建并反标准化，显著平滑测试链路的 NMSE/SNR 曲线并生成稳定的 CIR 输出
"""
import os
import math
import random
import numpy as np
from scipy.io import loadmat, savemat
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader
from tqdm import tqdm
import hdf5storage, h5py
from scipy.fft import dct

# =========================
# 1) 随机性与设备
def set_seed(seed=2025):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
set_seed(2025)

# =========================
# 2) 载入与预处理
# def load_cir_mat(mat_path, key='CIR'):
#     mat = loadmat(mat_path)
#     CIR = mat[key]  # 形状: (Nlink, Nrx, Ntx, Npath, Ntime)，复数
#     assert np.iscomplexobj(CIR), "CIR 应该是复数数组"
#     return CIR

def load_cir_mat(mat_path: str, cir_key: str = "CIR"):
    import numpy as np
    import hdf5storage, h5py

    # 只取指定变量名，避免把整个大文件一次性读进内存
    try:
        mat = hdf5storage.loadmat(mat_path, variable_names=[cir_key])
        if cir_key in mat:
            return np.asarray(mat[cir_key])
    except Exception as e:
        # 如果读取本身失败，继续列出键名辅助排错
        pass

    # 变量名不匹配或读取失败，列出可用顶层变量名
    with h5py.File(mat_path, "r") as f:
        keys = list(f.keys())
    raise KeyError(f"MAT 文件中找不到变量 '{cir_key}'。可用变量：{keys}")


def flatten_cir_to_sequences(CIR):
    """
    把 CIR [Nlink, Nrx, Ntx, Npath, Ntime] → 每条 link 生成一个时间序列：
    X_link: [Ntime, D]，其中 D = 2 * (Nrx*Ntx*Npath)
    复数拆成 [Re, Im] 级联。
    返回:
      X: list of np.ndarray with shape [T, D]
    """
    Nlink, Nrx, Ntx, Npath, Ntime = CIR.shape
    D = 2 * (Nrx * Ntx * Npath)
    X = []
    for i in range(Nlink):
        # reshape 成 [Nrx*Ntx*Npath, Ntime]
        mat_i = CIR[i].reshape(Nrx * Ntx * Npath, Ntime)
        re = np.real(mat_i).T  # [T, K]
        im = np.imag(mat_i).T
        x_seq = np.concatenate([re, im], axis=1)  # [T, 2K] = [T, D]
        X.append(x_seq.astype(np.float32))
    return X  # len = Nlink, 每项 [T, D]

# 标准化器：对D维向量做全局标准化（跨link+time）
class GlobalScaler:
    def __init__(self):
        self.scaler = StandardScaler(with_mean=True, with_std=True)

    def fit(self, X_list):
        """只在训练集数据上拟合标准化器以避免信息泄漏."""
        # 拼接训练集里所有 link / 时间步
        X_all = np.concatenate(X_list, axis=0)  # [sum_T_over_links, D]
        self.scaler.fit(X_all)
        return self

    def transform(self, X_list):
        return [self.scaler.transform(x) for x in X_list]

    def inverse_transform(self, X_list):
        return [self.scaler.inverse_transform(x) for x in X_list]

# =========================
# 3) 压缩测量与 PCA 粗重建
def build_measurement_matrix(D, M):
    # 原先返回 dense Phi 的实现会炸内存；改为返回隐式投影器
    return FastJLProjector(D, M, seed=2025)

def project_measurements(X_list, Phi):
    # 兼容：Phi 可能是隐式投影器，也可能是 ndarray
    if hasattr(Phi, "project"):
        return Phi.project(X_list)
    else:
        return [x.dot(Phi.T) for x in X_list]


class FastJLProjector:
    """
    隐式随机投影：y = S * DCT( diag(sign) * x )
    - 不显式生成 MxD 的 Phi
    - 仅存两样东西：随机 ±1 向量 sign[D] 和 选取的行索引 rows[M]
    - 对任意 [T, D] 的 x，按特征维做 DCT 后，抽取 M 列作为测量
    """
    def __init__(self, D, M, seed=2025):
        rng = np.random.default_rng(seed)
        self.D = int(D)
        self.M = int(M)
        self.sign = rng.choice([1.0, -1.0], size=self.D).astype(np.float32)
        # 无放回采样 M 个频域索引；如需可换成分层或低频优先
        self.rows = rng.choice(self.D, size=self.M, replace=False)

    def project_one(self, xTD: np.ndarray) -> np.ndarray:
        # xTD: [T, D] (float32/float64)，按列(轴=1)投影
        # 1) 随机符号翻转
        z = xTD * self.sign  # 广播做逐列 ±1
        # 2) 沿特征维做 DCT-II，正交归一
        Zf = dct(z, type=2, norm='ortho', axis=1)
        # 3) 抽取 M 列作为测量
        y = Zf[:, self.rows]
        # 4) 能量补偿：抽取列后乘 √(D/M)，匹配原向量均方能量
        y = y * math.sqrt(self.D / self.M)
        return y.astype(np.float32, copy=False)

    def project(self, X_list):
        return [self.project_one(x) for x in X_list]


class PCABranch:
    """
    用 PCA 学到稳态低秩子空间。拟合在标准化后的 X 上。
    提供：
      - transform：得到主成分系数 z
      - coarse_recon：PCA 逆变换得到低秩粗重建 x̂_pca
    """
    def __init__(self, n_components=64, random_state=2025):
        self.pca = PCA(n_components=n_components, svd_solver='auto', random_state=random_state)

    def fit(self, X_list_std):
        X_all = np.concatenate(X_list_std, axis=0)
        self.pca.fit(X_all)
        return self

    def transform(self, X_list_std):
        return [self.pca.transform(x) for x in X_list_std]

    def coarse_recon(self, X_list_std):
        return [self.pca.inverse_transform(self.pca.transform(x)) for x in X_list_std]

# =========================
# 4) BiGRU 残差学习（seq2seq）
class ResidualBiGRU(nn.Module):
    """
    输入:
      - y_seq  ∈ R^{T×M}   压缩测量序列
      - xp_seq ∈ R^{T×D}   PCA 粗重建（标准化域）
    输出:
      - r_hat  ∈ R^{T×D}   对标准化残差的估计
    """
    def __init__(self, M, D, hidden=256, num_layers=2, dropout=0.1):
        super().__init__()
        # BiGRU 的输入包含测量 y 以及 PCA 粗重建 xp
        input_size = M + D
        self.gru = nn.GRU(
            input_size=input_size,
            hidden_size=hidden,
            num_layers=num_layers,
            batch_first=True,
            dropout=dropout,
            bidirectional=True
        )
        self.head = nn.Sequential(
            nn.LayerNorm(2*hidden),
            nn.Linear(2*hidden, 2*hidden),
            nn.GELU(),
            nn.Linear(2*hidden, D)
        )

    def forward(self, y_seq, xp_seq):   # [B, T, M], [B, T, D]
        # 拼接观测与 PCA 粗重建，让网络利用两者的互补信息
        inp = torch.cat([y_seq, xp_seq], dim=-1)
        h, _ = self.gru(inp)  # [B, T, 2H]
        r = self.head(h)        # [B, T, D]
        return r

# 数据集：改为“滑动窗口小片段”，每个样本是一段长度为 win_len 的片段
class SlidingWindowDataset(Dataset):
    """
    从多条 link 的序列中抽取固定长度的滑窗片段：
      输入列表元素形状：
        Y_list[i]      : [T, M]
        Xstd_list[i]   : [T, D]
        Xpca_list[i]   : [T, D]
      返回单个样本：
        y_win  [win_len, M]
        xp_win [win_len, D] (PCA 粗重建)
        x_win  [win_len, D] (标准化真值)

    说明
    ----
    在某些旧版脚本里 ``LinkSequenceDataset`` 只接受三个位置参数
    ``(Y_list, Xstd_list, Xpca_list)``，不支持 ``win_len``/``stride`` 关键字。
    为了向后兼容，这里把窗口长度参数做成可选：当 ``win_len`` 为空或
    非正数时，会退化为“整段样本”模式，即每条链路仅产生一个样本。
    """

    def __init__(self, Y_list, Xstd_list, Xpca_list,
                 win_len=256, stride=128, dtype=np.float32, **_unused):
        assert len(Y_list) == len(Xstd_list) == len(Xpca_list)
        self.Y_list = Y_list
        self.Xstd_list = Xstd_list
        self.Xpca_list = Xpca_list
        self.dtype = dtype

        # ``win_len`` 或 ``stride``<=0 时退化为整段模式（兼容旧实现）。
        if win_len is None or win_len <= 0:
            self.win_len = None
            self.stride = None
        else:
            self.win_len = int(win_len)
            self.stride = int(stride if stride is not None else win_len)

        # 建立 (link_id, start) 索引。如果 win_len 为空，则每个 link 只保留整段。
        self.index = []
        for lid, (Y, X, XP) in enumerate(zip(Y_list, Xstd_list, Xpca_list)):
            T = Y.shape[0]
            if self.win_len is None:
                # 整段样本，直接把起点记为 0
                if T == 0:
                    continue
                self.index.append((lid, 0))
                continue

            if T < self.win_len:
                # 太短则跳过该 link（如需保留可改为 padding）
                continue
            s = 0
            while s + self.win_len <= T:
                self.index.append((lid, s))
                s += self.stride
            # 尾部对齐：确保覆盖到序列末尾
            if s < T:
                self.index.append((lid, T - self.win_len))

    def __len__(self):
        return len(self.index)

    def __getitem__(self, idx):
        lid, s = self.index[idx]
        if self.win_len is None:
            y = self.Y_list[lid].astype(self.dtype, copy=False)
            x = self.Xstd_list[lid].astype(self.dtype, copy=False)
            xp = self.Xpca_list[lid].astype(self.dtype, copy=False)
        else:
            e = s + self.win_len
            y  = self.Y_list[lid][s:e].astype(self.dtype, copy=False)
            x  = self.Xstd_list[lid][s:e].astype(self.dtype, copy=False)
            xp = self.Xpca_list[lid][s:e].astype(self.dtype, copy=False)
        # 返回测量、PCA 粗重建以及真值，方便外部计算残差/损失
        return (
            torch.from_numpy(y),
            torch.from_numpy(xp),
            torch.from_numpy(x)
        )


# ``LinkSequenceDataset`` 是旧版名称，保留别名以兼容外部引用。
LinkSequenceDataset = SlidingWindowDataset

class LinkFullSequenceDataset(Dataset):
    """
    测试/验证用：一个样本 = 一条 link 的完整序列，返回 (y_seq, xp_seq, x_seq)
    """
    def __init__(self, Y_list, Xstd_list, Xpca_list, dtype=np.float32):
        assert len(Y_list) == len(Xstd_list) == len(Xpca_list)
        self.Y_list = Y_list
        self.Xstd_list = Xstd_list
        self.Xpca_list = Xpca_list
        self.dtype = dtype

    def __len__(self):
        return len(self.Y_list)

    def __getitem__(self, idx):
        y  = self.Y_list[idx].astype(self.dtype, copy=False)
        x  = self.Xstd_list[idx].astype(self.dtype, copy=False)
        xp = self.Xpca_list[idx].astype(self.dtype, copy=False)
        return (
            torch.from_numpy(y),
            torch.from_numpy(xp),
            torch.from_numpy(x)
        )


def reconstruct_sequence_with_overlap(model, y_seq, xp_seq,
                                      win_len=None, stride=None):
    """在推理阶段使用与训练一致的滑窗策略，缓解长序列漂移导致的周期性误差。

    参数
    ----
    model : nn.Module
        训练好的 ``ResidualBiGRU``。
    y_seq : np.ndarray, shape [T, M]
        单条链路的压缩观测序列。
    xp_seq : np.ndarray, shape [T, D]
        与之对应的 PCA 粗重建（已标准化）。
    win_len : int or None
        训练时的窗口长度；若为空或不合法，则退化为整段一次性推理。
    stride : int or None
        滑窗步长；缺省时默认为 ``win_len``（即不重叠）。

    返回
    ----
    np.ndarray, shape [T, D]
        通过重叠-平均得到的标准化域重建结果。

    说明
    ----
    训练阶段的 DataLoader 每个样本都是长度 ``win_len`` 的片段，BiGRU 的隐状态
    在窗口边界被频繁重置。如果推理时一次性喂入远长于 ``win_len`` 的序列，
    隐状态会进入训练中未出现的分布，造成图中所示的周期性起伏。为保持训练/
    测试分布一致，这里将测试序列也拆成滑窗片段，模型在各片段内独立前向，
    然后使用 Hann 加权做重叠平均，以平滑窗口边界的拼接误差。
    """

    y_seq = np.asarray(y_seq, dtype=np.float32)
    xp_seq = np.asarray(xp_seq, dtype=np.float32)
    T = y_seq.shape[0]

    if win_len is None or win_len <= 0 or win_len >= T:
        # 训练中没有使用滑窗（或窗口比序列还长），直接整段推理
        with torch.no_grad():
            y_tensor = torch.from_numpy(y_seq).unsqueeze(0).to(device)
            xp_tensor = torch.from_numpy(xp_seq).unsqueeze(0).to(device)
            r_hat = model(y_tensor, xp_tensor).cpu().numpy()[0]
        return xp_seq + r_hat

    win_len = int(win_len)
    stride = int(stride) if (stride is not None and stride > 0) else win_len

    # 为重叠区域准备 Hann 权重，避免窗口边缘跳变
    weight_1d = np.hanning(win_len).astype(np.float32)
    if not np.any(weight_1d > 0):
        weight_1d = np.ones(win_len, dtype=np.float32)
    else:
        weight_1d = np.clip(weight_1d, 1e-3, None)
    weight_1d = weight_1d.reshape(win_len, 1)

    accum = np.zeros_like(xp_seq, dtype=np.float32)
    weight_sum = np.zeros((T, 1), dtype=np.float32)

    starts = []
    s = 0
    while s + win_len <= T:
        starts.append(s)
        s += stride
    if not starts:
        starts = [max(0, T - win_len)]
    elif starts[-1] != T - win_len:
        starts.append(max(0, T - win_len))

    for s in starts:
        e = s + win_len
        y_win = torch.from_numpy(y_seq[s:e]).unsqueeze(0).to(device)
        xp_win = torch.from_numpy(xp_seq[s:e]).unsqueeze(0).to(device)
        with torch.no_grad():
            r_win = model(y_win, xp_win).cpu().numpy()[0]
        x_win_hat = xp_seq[s:e] + r_win
        accum[s:e] += x_win_hat * weight_1d
        weight_sum[s:e] += weight_1d

    # 归一化重叠区域的权重和
    np.maximum(weight_sum, 1e-6, out=weight_sum)
    x_std_hat = accum / weight_sum
    return x_std_hat.astype(np.float32, copy=False)


# 训练循环
def train_bigru(model, loader, epochs=100, lr=1e-3, clip=1.0):
    """训练残差预测用的 BiGRU 网络。

    Parameters
    ----------
    model : nn.Module
        ``ResidualBiGRU`` 或兼容的模型，要求 `forward` 接受 ``(y_seq, xp_seq)``。
    loader : DataLoader
        迭代器会返回 ``(y, x_pca, x_std)`` 三元组，便于在循环内构造残差标签
        ``r_gt = x_std - x_pca``。
    epochs : int, optional
        训练轮数。
    lr : float, optional
        AdamW 的初始学习率。
    clip : float or None, optional
        梯度裁剪阈值；设置为 ``None`` 时跳过裁剪。

    Notes
    -----
    旧版本只把压缩观测 ``y`` 送入网络，模型无法利用 PCA 粗重建的先验信息。
    当前实现沿用原有接口，但会把 ``y`` 与 ``x_pca`` 一起输入模型，因此
    DataLoader 必须提供三元组批次。这样 BiGRU 在学习残差时能够感知 PCA
    子空间的粗重建，为最终的复原提供额外上下文。
    """
    model = model.to(device)
    opt = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=1e-4)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=epochs)
    loss_fn = nn.SmoothL1Loss()  # 对残差更稳健
    best = {'loss': float('inf'), 'state': None}

    for ep in range(1, epochs+1):
        model.train()
        running = 0.0
        for y, xp, x in loader:
            y = y.to(device)         # [B, T, M]
            xp = xp.to(device)       # [B, T, D]
            x = x.to(device)         # [B, T, D]
            opt.zero_grad()
            r_gt = x - xp            # 真实残差（标准化域）
            r_hat = model(y, xp)
            loss = loss_fn(r_hat, r_gt)
            loss.backward()
            if clip is not None:
                nn.utils.clip_grad_norm_(model.parameters(), clip)
            opt.step()
            running += loss.item() * y.size(0)
        scheduler.step()
        ep_loss = running / len(loader.dataset)
        print(f"[Epoch {ep:03d}] train loss: {ep_loss:.6f}")
        if ep_loss < best['loss']:
            best['loss'] = ep_loss
            best['state'] = {k: v.cpu().clone() for k, v in model.state_dict().items()}

    if best['state'] is not None:
        model.load_state_dict(best['state'])
    return model

# =========================
# 5) 评估指标与重建流程
# =========================
def nmse(true, pred, eps=1e-12):
    # true/pred: [T, D]
    num = np.sum((true - pred) ** 2)
    den = np.sum(true ** 2) + eps
    return num / den

def recon_snr(true, pred, eps=1e-12):
    # 10*log10( ||x||^2 / ||x - x̂||^2 )
    num = np.sum(true ** 2)
    den = np.sum((true - pred) ** 2) + eps
    return 10.0 * np.log10(num / den)

def run_pipeline(
    mat_path='C:/Users/admin/Desktop/gen_dd_channel/result_Racian/CIR.mat',
    cir_key='CIR',
    out_dir='./outputs4',
    pca_components=64,
    compress_ratio=0.25,
    win_len=256, win_stride=128,   # ★ 新增
    batch_size=2, epochs=50, lr=1e-3
):

    global t_split
    t_split = None  # 默认 None，单链路时再赋值
    os.makedirs(out_dir, exist_ok=True)

    # ---- 载入 CIR 并展平 ----
    CIR = load_cir_mat(mat_path, cir_key)              # complex
    X_list = flatten_cir_to_sequences(CIR)             # list of [T, D]
    T0, D = X_list[0].shape
    print(f"#Links={len(X_list)}, T={T0}, D={D}")

    # ---- 划分训练/验证（先划分，避免预处理时泄漏） ----
    Nlink = len(X_list)

    if Nlink >= 2:
        # 多链路（保留：即便你现在是单链路，将来也能复用）
        link_idx = np.arange(Nlink)
        tr_idx, te_idx = train_test_split(link_idx, test_size=0.3, random_state=2025)
    else:
        # 单链路：按时间划分（70% 训练，30% 测试）
        T = X_list[0].shape[0]
        if T < 2:
            raise ValueError("单链路序列长度至少需要 2 个时间步以划分训练/测试。")
        t_split = int(T * 0.7)
        t_split = max(1, min(t_split, T - 1))
        tr_idx = [0]
        te_idx = [0]

    # ---- 标准化：仅用训练子集拟合，随后变换整库 ----
    X_train_for_fit = [X_list[i] for i in tr_idx] if Nlink >= 2 else [X_list[0][:t_split]]
    scaler = GlobalScaler().fit(X_train_for_fit)
    X_std = scaler.transform(X_list)

    # ---- PCA：同样只用训练子集拟合 ----
    X_train_std_for_pca = [X_std[i] for i in tr_idx] if Nlink >= 2 else [X_std[0][:t_split]]
    pca = PCABranch(n_components=min(pca_components, D//2)).fit(X_train_std_for_pca)
    X_pca_std = pca.coarse_recon(X_std)

    # ---- 测量矩阵与压缩测量 ----
    M = max(8, int(compress_ratio * D))
    Phi = build_measurement_matrix(D, M)               # 隐式随机投影
    Y_list = project_measurements(X_std, Phi)          # list of [T, M]

    if Nlink >= 2:
        # 训练：滑窗小片段
        train_ds = SlidingWindowDataset(
            [Y_list[i] for i in tr_idx],
            [X_std[i] for i in tr_idx],
            [X_pca_std[i] for i in tr_idx],
            win_len=win_len, stride=win_stride
        )

        # 测试：整段
        test_ds = LinkFullSequenceDataset(
            [Y_list[i] for i in te_idx],
            [X_std[i] for i in te_idx],
            [X_pca_std[i] for i in te_idx]
        )
    else:
        # 训练：滑窗小片段（只用前 70%）
        train_ds = SlidingWindowDataset(
            [Y_list[0][:t_split]],
            [X_std[0][:t_split]],
            [X_pca_std[0][:t_split]],
            win_len=win_len, stride=win_stride
        )
        # 若训练片段太短导致没有样本，兜底缩小窗口
        if len(train_ds) == 0:
            adj_win = min(win_len, max(8, t_split))
            train_ds = SlidingWindowDataset(
                [Y_list[0][:t_split]],
                [X_std[0][:t_split]],
                [X_pca_std[0][:t_split]],
                win_len=adj_win, stride=adj_win
            )

        # 测试：整段（用后 30%）
        test_ds = LinkFullSequenceDataset(
            [Y_list[0][t_split:]],
            [X_std[0][t_split:]],
            [X_pca_std[0][t_split:]]
        )

    train_loader = DataLoader(
        train_ds, batch_size=batch_size, shuffle=True,
        drop_last=False, num_workers=0, pin_memory=False
    )
    test_loader = DataLoader(
        test_ds, batch_size=1, shuffle=False,
        drop_last=False, num_workers=0, pin_memory=False
    )

    # ---- 训练 BiGRU 残差网络 ----
    model = ResidualBiGRU(M=M, D=D, hidden=256, num_layers=2, dropout=0.1)
    model = train_bigru(model, train_loader, epochs=epochs, lr=lr)

    # ---- 推理与重建（回到原尺度）----
    model.eval()
    nmse_list, snr_list = [], []
    recon_all = []   # 存放每个测试 link 的重建（原尺度）

    eval_win = getattr(train_ds, 'win_len', win_len)
    eval_stride = getattr(train_ds, 'stride', win_stride)

    for (y, xp, x_std_gt) in tqdm(test_loader, desc="Evaluate"):
        # DataLoader 默认返回 CPU Tensor，直接转换为 numpy
        y_np = y.numpy()[0]
        xp_np = xp.numpy()[0]
        x_std_gt_np = x_std_gt.numpy()[0]

        x_std_hat = reconstruct_sequence_with_overlap(
            model, y_np, xp_np, win_len=eval_win, stride=eval_stride
        )

        # 反标准化回原值域
        x_hat = scaler.scaler.inverse_transform(x_std_hat)
        x_gt  = scaler.scaler.inverse_transform(x_std_gt_np)

        # 评估
        nmse_val = nmse(x_gt, x_hat)
        snr_val  = recon_snr(x_gt, x_hat)
        nmse_list.append(nmse_val)
        snr_list.append(snr_val)
        recon_all.append(x_hat.astype(np.float32))

    print(f"[Test]  Avg NMSE = {np.mean(nmse_list):.6f}")
    print(f"[Test]  Avg SNR  = {np.mean(snr_list):.2f} dB")

    # ---- 把向量重构回原复数 CIR 格式 ----
    # 原始每个时间步的向量 x(t) = [Re(vec); Im(vec)], vec 长度 K = Nrx*Ntx*Npath
    Nlink, Nrx, Ntx, Npath, Ntime = CIR.shape
    K = Nrx * Ntx * Npath
    # 只把“测试集的 link”重建回复数 CIR 形状，保存 npy
    recon_cir_links = []
    for x_hat in recon_all:  # [T, D]
        re = x_hat[:, :K]
        im = x_hat[:, K:]
        vec = re + 1j * im           # [T, K]
        ci = vec.T.reshape(Nrx, Ntx, Npath, -1)  # [Nrx, Ntx, Npath, T]
        recon_cir_links.append(ci)

    np.save(os.path.join(out_dir, 'recon_cir.npy'), recon_cir_links, allow_pickle=True)
    # savemat(os.path.join(out_dir, 'recon_cir.mat'), {'recon_cir_links': np.array(recon_cir_links, dtype=object)})

    # 同时返回关键对象，方便继续实验
    return {
        'Phi': Phi, 'PCA': pca, 'Scaler': scaler,
        'model': model, 'metrics': {'NMSE': nmse_list, 'SNR': snr_list},
        'recon_links': recon_cir_links,
        't_split': t_split
    }

if __name__ == "__main__":
    result = run_pipeline(
        mat_path='C:/Users/admin/Desktop/gen_dd_channel/result_Racian/CIR.mat',
        cir_key='CIR',
        out_dir='./outputs4',
        pca_components=64,  # 可调：PCA 低秩维度
        compress_ratio=0.25,  # 可调：压缩比 M/D
        batch_size=2,
        epochs=50,
        lr=1e-3
    )

# =================================
# 5).可视化
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat

# ---- 配置：路径 & 选择的链路/Rx/Tx ----
MAT_PATH = r"C:\Users\admin\Desktop\gen_dd_channel\result_Racian\CIR.mat"
MAT_KEY = "CIR"
# RECON_NPY = r"D:\PythonProject\outputs2\recon_cir.npy"  # 重建输出
RECON_NPY = os.path.abspath(os.path.join("outputs4", "recon_cir.npy"))
OUTDIR = r"D:\PythonProject\viz_out4"
LINK_IN_MAT = 0
LINK_IN_RECON = 0
IRX, ITX = 1, 1  # 画哪个 Rx/Tx
TOPK = 16  # 只取能量最高的 K 个 tap，可设 None 画全部

os.makedirs(OUTDIR, exist_ok=True)


def _ensure_5d(C):
    if C.ndim == 4: C = C[np.newaxis, ...]
    return C


def _pick_stream(C, irx=0, itx=0):  # [Nrx,Ntx,Npath,T] -> [Npath,T]
    return C[irx, itx, :, :]


def _nmse_t(true, pred, eps=1e-12):
    num = np.sum(np.abs(true - pred) ** 2, axis=0)
    den = np.sum(np.abs(true) ** 2, axis=0) + eps
    return num / den


def _snr_t(true, pred, eps=1e-12):
    num = np.sum(np.abs(true) ** 2, axis=0) + eps
    den = np.sum(np.abs(true - pred) ** 2, axis=0) + eps
    return 10 * np.log10(num / den)


# ---- 载入真值与重建（全程只用 4 维）----
mat_arr = load_cir_mat(MAT_PATH, MAT_KEY)
# 选出 4 维的真值 C_true
if mat_arr.ndim == 5:
    # 多链路：选定 LINK_IN_MAT
    if LINK_IN_MAT >= mat_arr.shape[0]:
        print(f"[warn] LINK_IN_MAT={LINK_IN_MAT} 越界，只有 {mat_arr.shape[0]} 条链路，改为 0")
        LINK_IN_MAT = 0
    C_true = mat_arr[LINK_IN_MAT]  # -> [Nrx,Ntx,Npath,T_true]
elif mat_arr.ndim == 4:
    # 单链路：本身就是 4 维
    C_true = mat_arr
else:
    raise ValueError(f"Unexpected CIR ndim={mat_arr.ndim}, expected 4 or 5.")

# 重建：列表里的 4 维
recon_list = np.load(RECON_NPY, allow_pickle=True)  # list of [Nrx,Ntx,Npath,T_hat]
if LINK_IN_RECON >= len(recon_list):
    print(f"[warn] LINK_IN_RECON={LINK_IN_RECON} 越界，只有 {len(recon_list)}，改为 0")
    LINK_IN_RECON = 0
C_hat = recon_list[LINK_IN_RECON]  # [Nrx,Ntx,Npath,T_hat]
# 与重建长度对齐：只取真值末尾 T_hat 帧
T_true = C_true.shape[-1]
T_hat = C_hat.shape[-1]
if T_hat <= T_true:
    C_true = C_true[..., T_true - T_hat:]
else:
    raise ValueError(f"T_hat({T_hat}) > T_true({T_true})，请检查重建长度或数据。")

assert C_true.shape == C_hat.shape, f"shape mismatch: {C_true.shape} vs {C_hat.shape}"

H_true = _pick_stream(C_true, IRX, ITX)  # [Npath,T]
H_hat = _pick_stream(C_hat, IRX, ITX)
Npath, T = H_true.shape

# 只画能量最高的 K 个 tap（可选）
if TOPK is not None and 1 <= TOPK < Npath:
    pwr = np.mean(np.abs(H_true) ** 2, axis=1)
    keep = np.argsort(pwr)[::-1][:TOPK]
    H_true = H_true[keep, :]
    H_hat = H_hat[keep, :]
    Npath = TOPK

# ---- 图1：多径叠加的幅度 / 相位 ----
H_sum_true = np.sum(H_true, axis=0)  # [T]
H_sum_hat = np.sum(H_hat, axis=0)
t = np.arange(T)

plt.figure(figsize=(9, 4))
plt.plot(t, 20 * np.log10(np.abs(H_sum_true) + 1e-12), label='True')
plt.plot(t, 20 * np.log10(np.abs(H_sum_hat) + 1e-12), '--', label='Recon')
plt.xlabel('time index');
plt.ylabel('|∑h| (dB)');
plt.legend()
plt.title('Aggregate amplitude (multipath sum)')
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "sum_amp_time.png"))
plt.close()

plt.figure(figsize=(9, 4))
plt.plot(t, np.unwrap(np.angle(H_sum_true)), label='True')
plt.plot(t, np.unwrap(np.angle(H_sum_hat)), '--', label='Recon')
plt.xlabel('time index');
plt.ylabel('phase (rad)');
plt.legend()
plt.title('Aggregate phase (multipath sum)')
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "sum_phase_time.png"))
plt.close()

# ---- 图2：最强 tap 的幅度/相位随时间 ----
tap_sel = int(np.argmax(np.mean(np.abs(H_true) ** 2, axis=1)))

plt.figure(figsize=(9, 4))
plt.plot(t, 20 * np.log10(np.abs(H_true[tap_sel, :]) + 1e-12), label='True')
plt.plot(t, 20 * np.log10(np.abs(H_hat[tap_sel, :]) + 1e-12), '--', label='Recon')
plt.xlabel('time index');
plt.ylabel(f'|h| @ tap {tap_sel} (dB)');
plt.legend()
plt.title('Amplitude vs time (strongest tap)')
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "tap_amp_time.png"));
plt.close()

plt.figure(figsize=(9, 4))
plt.plot(t, np.unwrap(np.angle(H_true[tap_sel, :])), label='True')
plt.plot(t, np.unwrap(np.angle(H_hat[tap_sel, :])), '--', label='Recon')
plt.xlabel('time index');
plt.ylabel('phase (rad)');
plt.legend()
plt.title('Phase vs time (strongest tap)')
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "tap_phase_time.png"));
plt.close()

# ---- 图4：逐时刻 NMSE / SNR ----
nmse_curve = _nmse_t(H_true, H_hat)
snr_curve = _snr_t(H_true, H_hat)

plt.figure(figsize=(9, 4))
plt.plot(t, nmse_curve)
plt.xlabel('time index');
plt.ylabel('NMSE');
plt.title('NMSE over time')
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "nmse_over_time.png"));
plt.close()

plt.figure(figsize=(9, 4))
plt.plot(t, snr_curve)
plt.xlabel('time index');
plt.ylabel('SNR (dB)');
plt.title('Reconstruction SNR over time')
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "snr_over_time.png"));
plt.close()

print(f"[OK] figures saved to: {OUTDIR}")
