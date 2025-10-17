% 基于D 波段双向测量数据生成无线信道实现的函数
% 极化域被省略，因为测量数据中仅包含 V-V 极化。
% 
% 输入： D 波段双向测量数据，存放于目录 data_030123_double_directional_paths。
% 可选项：
% 用于传输函数生成的信道模型类型：
%'snapshot'：每条链路仅取一个时间样本 →静态；每条路径分配随机相位，幅度取自测量数据。
% 'virtualMotion'：引入多普勒频移和时间衰落，传播参数保持静态；每条路径随机相位，幅度取自测量数据，多普勒频率由到达角 (AoA) 和速度向量计算得到。
% 天线/波束形状：
% 'single3GPP'：单天线单元，方向图形状由 3GPP 定义，可调半功率波束宽度 (HPBW) 等。
% 'URA'：均匀矩形阵列，单元全向。
% 'UCA'：均匀圆阵列，单元全向。
% 
% 输出：
% 'CIR'：信道冲激响应，维度为 Nlink × NRx × NTx × Npath × Ntime
% 'CFR'：信道频率响应，维度为 Nlink × NRx × NTx × Nfreq × Ntime
% 可选的可视化图：天线、CIR 和 CFR 图示
clear; format compact;

%% Load data %%%%%%%%%%%%%%%%%%%
% %% Load data %%%%%%%%%%%%%%%%%%% 
% 载入测量数据
% 设置要使用的场景，可选项包括：
% 'EntranceHall' (入口大厅) 
% 'Suburban' (郊区)
% 'Residential' (居民区)
% 'CityCenter' (市中心)
scenario = 'Residential'; 

% 调用自定义函数 gen_pathDD() 读取并转换提供的D波段双向测量数据，返回一个包含所有链路的结构体数组
pathData = gen_pathDD(scenario);    % 1×Links 的结构体，包含多径参数


visualizeAntenna = 0;       % 是否绘制天线和波束相关图形（0=不画，1=画）
visualizeChModel = 1;       % 是否绘制生成的信道模型结果图（0=不画，1=画）
visualizeLayout  = 1;       % 是否绘制Tx/Rx和散射体位置布局图（0=不画，1=画）

% 选择信道模型类型：
% 'snapshot'：仅生成一个时间快照（静态信道）
% 'virtualMotion'：包含多普勒频移和时间衰落，生成时间变化信道
% 'continuous'：插值链路参数，得到连续时间变化信道（尚未实现）
par.modelType = 'virtualMotion';

% 设置发射机和接收机天线阵列类型：
% 'single3GPP'：单天线单元，3GPP标准方向图
% 'URA'：均匀矩形阵列
% 'UCA'：均匀圆阵列
par.tx.antennaType = 'URA';  % 发射端天线类型
par.rx.antennaType = 'UCA';  % 接收端天线类型

par.fc = 140E9;             % 仿真中心频率 [Hz]
par.BW = 1E9;               % 信道带宽 [Hz]

% 指定要生成的链路索引（即选择哪些Tx-Rx链路）
% 可根据测量数据文件中提供的链路数量进行选择
% par.linkSet = [1 40 ]; 
par.linkSet = [1];     % 当前选择链路索引

SD = 2;                     % 采样密度：每半波长的采样点数（过采样因子）
Nlink = length(par.linkSet); % 要建模的链路数量
SoL = 299792458;            % 光速 [m/s]
lambda = SoL/par.fc;        % 波长 [m]

% 防止出现NaN值，设定距离下限
par.d_min     = max(0.01, 0.1*lambda);   % >=1 cm 或 λ/10
par.amp_clip  = 1e3;                     % 幅度上限（60 dB）
par.amp_floor = 1e-3;                    % 幅度下限（-60 dB）

% 根据选择的模型类型设置时间采样与接收机运动参数
switch par.modelType
    case 'snapshot'      % 单时间快照（静态信道）
        par.Ntime = 1;   % 每条链路仅生成1个时间样本
    case 'virtualMotion' % 虚拟运动（包含多普勒）
        par.Ntime = 5000;   % 每条链路的时间样本数，决定生成多少个时间点，相当于仿真多长的时间序列
        par.speed = 10;     % 接收机运动速度 [m/s]决定多普勒频移的幅度
        par.DoTaz = 45;     % 运动方向的方位角 [deg]决定运动方向，影响多普勒是正偏移还是负偏移
        par.DoTel = -5;     % 运动方向的俯仰角 [deg]决定运动方向，影响多普勒是正偏移还是负偏移
    case 'continuous'   % 连续时间信道（尚未实现）
        % 这里预留了插值生成动态信道的接口
    otherwise
        error('Choose proper par.modelType parameter value!')
end

par  % 显示当前设置的所有参数，便于检查


%% Tx波束/天线阵列定义 %%%%%%%%%%%%%%%%%%%%
% 说明：在 ULA / UCA 两类阵列里默认各阵元为全向（omni）元素
switch par.tx.antennaType
    case 'single3GPP'
        Nant = 1;           % 天线元素数量（单元天线，必须为 1）

        % 3GPP 功率方向图（TR 38.901 表 7.3-1）
        az = -179:180;      % 采样的方位角网格 [deg]（步进 1 度，不要改）
        el = -90:90;        % 采样的俯仰角网格 [deg]（步进 1 度，不要改）
        theta_3dB = 20;     % 竖直方向（俯仰）半功率波束宽度 HPBW [deg]
        phi_3dB   = 10;     % 水平方向（方位）半功率波束宽度 HPBW [deg]
        Amax = 30;          % 最大旁瓣/截止衰减 [dB]
        G_Emax = 8;         % 单元最大方向性增益 [dBi]

        % V cut（俯仰方向）衰减曲线：min(12*(el/θ3dB)^2, Amax)
        Aver_dB = -min(12*(el/theta_3dB).^2, Amax);
        % H cut（方位方向）衰减曲线：min(12*(az/φ3dB)^2, Amax)
        Ahor_dB = -min(12*(az/phi_3dB).^2, Amax);

        % 组合二维方向图矩阵（Nelev x Nazim），并限制最大衰减 Amax
        % A_dB = G_Emax - min( -(A_H + A_V), Amax )
        A_dB = G_Emax - min( -(repmat(Ahor_dB, length(el), 1) ...
                             + repmat(Aver_dB', 1, length(az))), Amax);

        % 波束指向设置：'max' 表示自适应对准最大功率入射方向；
        % 也可固定成 [azi, ele] 度数
%         orientA = [23, -15];  % 若想手动指定指向角，取消注释并设置
        orientA = 'max';

        % 如果使用固定指向角，则把方向图按给定[俯仰,方位]做整度平移
        if ~strcmpi(orientA,'max')
            A_dB = circshift(A_dB, [orientA(2) orientA(1)]);
        end

        % 可视化天线方向图
        if visualizeAntenna
            figure(1)
            subplot(3,1,1)
            plot(az, G_Emax + Ahor_dB), grid on
            xlabel('azimuth [deg]'), ylabel('radiation power pattern [dB]')
            title(sprintf('Tx Single 3GPP beam, horizontal cut, HPBW = %.1f deg', phi_3dB))

            subplot(3,1,2)
            plot(el, G_Emax + Aver_dB), grid on
            xlabel('elevation [deg]'), ylabel('radiation power pattern [dB]')
            title(sprintf('Tx Single 3GPP beam, vertical cut HPBW = %.1f deg', theta_3dB))

            subplot(3,1,3)
            surf(az, el, A_dB, 'EdgeColor', 'flat'), axis equal, view(2)
            xlabel('azimuth [deg]'), ylabel('elevation [deg]'), zlabel('gain [dB]')
            colorbar
        end

        % 注意：3GPP 下方位范围为 -180:180，而测量数据使用 0:360，因此做一次平移
        A_dB = circshift(A_dB, [0 -180]);
        az = 1:360;              % 重定义 az 索引与上面的平移保持一致
        tx.A_dB = A_dB;          % 保存方向图（dB）
        tx.orientA = orientA;    % 保存指向设置
        tx.acoord = [0 0 0];     % 单元天线时仅用占位坐标

    case 'URA'
        % URA（Uniform Rectangular Array，均匀矩形阵）
        Nrow = 4;               % 阵列的行数（竖直方向）
        Ncol = 8;               % 阵列的列数（水平方向）
        Nant = Nrow * Ncol;     % 总阵元数
        rowSpacing = 0.5;       % 行间距（以波长为单位）
        colSpacing = 0.5;       % 列间距（以波长为单位）

        % 波束指向：'max' 表示对准最大功率路径；也可固定为具体 [azi, ele]
%         orientA = [23, -15];
        orientA = 'max';

        % 生成阵元三维坐标（单位：米）
        xa = zeros(1, Nant);                % x 坐标：平面阵，取 0
        ytmp = 0:colSpacing:colSpacing*(Ncol-1);
        ytmp = (ytmp - max(ytmp)/2);        % 水平坐标居中
        ztmp = 0:rowSpacing:rowSpacing*(Nrow-1);
        ztmp = (ztmp - max(ztmp)/2);        % 垂直坐标居中
        [ya, za] = meshgrid(ytmp, ztmp);
        ya = ya(:)' * lambda;               % 转为米（乘以波长）
        za = za(:)' * lambda;               % 转为米（乘以波长）
        acoord = [xa', ya', za'];           % 阵元坐标矩阵 [Nant x 3]

        % 如采用固定指向角，则把阵列刚体旋转到该指向
        if ~strcmpi(orientA,'max')
            rotEle = roty(orientA(2));      % 绕 y 轴的俯仰旋转矩阵
            rotAzi = rotz(orientA(1));      % 绕 z 轴的方位旋转矩阵
            tmp = rotEle * acoord';         % 先做俯仰旋转
            tmp = rotAzi * tmp;             % 再做方位旋转
            acoord = tmp';
        end
        tx.acoord = acoord;     % 保存阵元坐标
        tx.orientA = orientA;   % 保存指向设置

    case 'UCA'
        % UCA（Uniform Circular Array，均匀圆阵）
        % 由于仅用到方位域数据，不必定义波束正前方向（boresight）
        aper = 8;                           % 圆阵口径（直径，单位：波长）
        Nant = ceil(aper * pi * 2);         % 以 <= λ/2 的间距放置阵元,计算元素数
        [xa, ya] = pol2cart(linspace(0, 2*pi, Nant+1), aper/2 * lambda);
        acoord = [xa', ya', zeros(length(xa), 1)];  % 阵元坐标（位于 xy 平面）
        acoord(end,:) = [];                 % 移除最后一个重复点
        tx.acoord = acoord;                 % 保存阵元坐标
        tx.orientA = 'uca';                 % 标记为圆阵

    otherwise
        error('Check the par.antennaType parameter value!') % 非法类型报错
end

tx.Nant = Nant;  % 记录 Tx 阵元数

% 可视化阵列综合方向图（仅对 URA/UCA 可视）
if visualizeAntenna & sum(strcmpi(par.tx.antennaType, {'URA','UCA'}))
    figure(1)
    subplot(3,2,5)

    % 在一组 az/el 网格上扫描波束（计算阵列因子）
    [azq, elq] = meshgrid([-120:1:120], [-60:1:60]);  % 扫描角度范围
    % 根据扫描角生成对应的波矢（单位：1/m）
    betaq = [cosd(azq(:)).*cosd(elq(:)), ...
             sind(azq(:)).*cosd(elq(:)), ...
             sind(elq(:))] * 2*pi/lambda;

    % 设定希望的波束指向（用于生成数字波束权向量）
    angAzSt = 0;                         % 目标指向：方位角 0°
    angElSt = 0;                         % 目标指向：俯仰角 0°
    % 指向 (angAzSt, angElSt) 的波矢
    gamma = -[cosd(angAzSt).*cosd(angElSt), ...
              sind(angAzSt).*cosd(angElSt), ...
              sind(angElSt)] * 2*pi/lambda;
    % 对应的理想加权（steering vector）
    vecSt = exp(1j * gamma * acoord');   % [1 x Nant] 复权值

    % 计算在各个扫描角上的阵列因子（含期望指向的相位对齐）
    af = sum(exp(1j * (betaq * acoord' + repmat(angle(vecSt), length(betaq), 1))), 2);
    afMat = reshape(af, size(azq,1), size(azq,2));  % 变成 (elev x azi) 网格

    % 绘制二维热力图（阵列因子幅度）
    surf(azq, elq, 20*log10(abs(afMat)), 'EdgeColor', 'flat'),
    axis equal, view(2), colorbar
    xlabel('azimuth [deg]'), ylabel('elevation [deg]'), zlabel('array factor [dB]')

    % 水平切片（el=0°）的阵列因子曲线
    subplot(3,2,[1 2])
    ind = find(elq == 0);
    plot(azq(ind), 20*log10(abs(afMat(ind)))), grid on
    set(gca, 'XTick', [-135:45:135]), xlim([-135 135])
    xlabel('azimuth [deg]'), ylabel('array factor [dB]')
    title(sprintf('Tx Beam steered to az=%.1f,  el=%.1f, horizontal cut at 0 elevation', angAzSt, angElSt))

    % 竖直切片（az=0°）的阵列因子曲线
    subplot(3,2,[3 4])
    ind = find(azq == 0);
    plot(elq(ind), 20*log10(abs(afMat(ind)))), grid on
    xlabel('elevation [deg]'), ylabel('array factor [dB]')
    title(sprintf('Tx Beam steered to az=%.1f,  el=%.1f, horizontal cut at 0 elevation', angAzSt, angElSt))

    % 绘制阵列的几何位置（单位：mm）
    subplot(3,2,6)
    plot3(tx.acoord(:,1)*1E3, tx.acoord(:,2)*1E3, tx.acoord(:,3)*1E3, 'o')
    grid on, axis square equal
    xlabel('X [mm]'), ylabel('Y [mm]'), zlabel('Z [mm]')
    title(sprintf('Tx Array geometry, %d element %s', Nant, par.tx.antennaType))
end



%% 定义接收端波束/天线阵列 %%%%%%%%%%%%%%%%%%%% 
% 选择 URA/UCA时调用，假设每个天线单元都是理想全向天线
switch par.rx.antennaType
    case 'single3GPP'  % 情况1：单天线，使用3GPP定义的方向图
        Nant = 1;           % 天线元素数量 = 1
        % 载入3GPP标准中定义的辐射方向图功率模式（TR 38.901表7.3-1）
        az = -179:180;      % 方位角采样范围[-179°~180°]，步长1°
        el = -90:90;        % 俯仰角采样范围[-90°~90°]，步长1°
        theta_3dB = 20;     % 垂直方向半功率波束宽度 (HPBW) [°]
        phi_3dB = 10;       % 水平方向半功率波束宽度 (HPBW) [°]
        Amax = 30;          % 最大旁瓣衰减 [dB]
        G_Emax = 8;         % 天线最大方向增益 [dBi]

        % 计算垂直方向增益曲线：根据3GPP经验公式，抛物线收敛到Amax
        Aver_dB = -min(12*(el/theta_3dB).^2,Amax);  
        % 计算水平方向增益曲线
        Ahor_dB = -min(12*(az/phi_3dB).^2,Amax);    
        % 生成二维增益矩阵 (俯仰x方位)
        A_dB = G_Emax - min(-(repmat(Ahor_dB,length(el),1) ...
                           + repmat(Aver_dB',1,length(az))),Amax);
        
        % 设置波束指向：'max' 表示自动对准最大功率路径，也可固定角度
        orientA = 'max';    

        % 如果指定了固定角度，就对增益矩阵进行循环移位以实现旋转
        if ~strcmpi(orientA,'max')
            A_dB = circshift(A_dB,[orientA(2) orientA(1)]);
        end
        
        % 可视化天线方向图（水平切片、垂直切片、二维辐射图）
        if visualizeAntenna
            figure(3)
            subplot(3,1,1)
            plot(az,G_Emax+Ahor_dB), grid on
            xlabel('azimuth [deg]'),ylabel('radiation power pattern [dB]')
            title(sprintf('Rx 单天线水平方向图, HPBW = %.1f deg',phi_3dB))
            subplot(3,1,2)
            plot(el,G_Emax+Aver_dB), grid on
            xlabel('elevation [deg]'),ylabel('radiation power pattern [dB]')
            title(sprintf('Rx 单天线垂直方向图, HPBW = %.1f deg',theta_3dB))
            subplot(3,1,3)
            surf(az,el,A_dB, 'EdgeColor','flat'), axis equal, view(2)
            xlabel('azimuth [deg]'),ylabel('elevation [deg]'),zlabel('gain [dB]')
            colorbar
        end

        % 将方向图矩阵对齐到0~360°坐标（因为原始数据用0~360°）
        A_dB = circshift(A_dB,[0 -180]);
        az = 1:360;
        rx.A_dB = A_dB;      % 存储方向图
        rx.orientA = orientA;
        rx.acoord = [0 0 0]; % 单天线坐标（仅占位）

    case 'URA' % 接收端为均匀矩形阵列
        Nrow = 2;           % 阵列行数
        Ncol = 3;           % 阵列列数
        Nant = Nrow*Ncol;   % 总天线数
        rowSpacing = 0.5;   % 行间距 = 0.5波长
        colSpacing = 0.5;   % 列间距 = 0.5波长
        orientA = 'max';    % 波束指向最大路径方向（可改为固定角度）

        % 生成URA元素在三维空间的坐标
        xa = zeros(1,Nant);  % x坐标（假设放在xz平面上）
        ytmp = 0:colSpacing:colSpacing*(Ncol-1); 
        ytmp = (ytmp - max(ytmp)/2); % 居中
        ztmp = 0:rowSpacing:rowSpacing*(Nrow-1); 
        ztmp = (ztmp - max(ztmp)/2);
        [ya,za] = meshgrid(ytmp,ztmp);
        ya = ya(:)'*lambda;  % 转换为实际长度（乘波长）
        za = za(:)'*lambda;
        acoord = [xa', ya', za'];  % 所有天线的三维坐标
        
        % 如果指定了固定指向，对阵列坐标做旋转
        if ~strcmpi(orientA,'max')
            rotEle = roty(orientA(2));  % 绕y轴旋转
            rotAzi = rotz(orientA(1));  % 绕z轴旋转
            tmp = rotEle*acoord';       
            tmp = rotAzi*tmp;           
            acoord = tmp';
        end
        rx.acoord = acoord;  % 存储URA坐标
        rx.orientA = orientA;

    case 'UCA' % 情况3：接收端为均匀圆形阵列
        aper = 8;                % 圆阵直径（单位：波长）
        Nant = ceil(aper*pi*2);  % 确定元素数目（间距≤λ/2）
        [xa,ya] = pol2cart(linspace(0,2*pi,Nant+1),aper/2*lambda);
        acoord = [xa', ya', zeros(length(xa),1)]; % 每个元素坐标
        acoord(end,:) = [];      % 删除重复的最后一个点
        rx.acoord = acoord;
        rx.orientA = 'uca';
        
    otherwise
        error('接收端天线类型参数设置错误！')
end

rx.Nant = Nant;  % 存储接收端总天线数

% URA、UCA可视化，绘制阵列方向图
if visualizeAntenna & sum(strcmpi(par.rx.antennaType,{'URA','UCA'})) 
    figure(2)
    subplot(3,2,5)
    % 生成扫描角度网格
    [azq,elq] = meshgrid([-120:1:120],[-60:1:60]);   
    betaq = [cosd(azq(:)).*cosd(elq(:)),...
             sind(azq(:)).*cosd(elq(:)),...
             sind(elq(:))]*2*pi/lambda; % 计算波矢
    
    % 生成波束指向方向的导向矢量
    angAzSt = 0; angElSt = 0; % 扫描角度设为0°指向
    gamma = -[cosd(angAzSt).*cosd(angElSt), sind(angAzSt).*cosd(angElSt), sind(angElSt)]*2*pi/lambda;
    vecSt = exp(1j*gamma*acoord'); 
    af = sum(exp(1j*(betaq*acoord' + repmat(angle(vecSt),length(betaq),1)) ),2); 
    afMat = reshape(af,size(azq,1),size(azq,2));  % 重塑为二维矩阵
    surf(azq,elq,20*log10(abs(afMat)),'EdgeColor','flat'), 
    axis equal, view(2), colorbar
    xlabel('azimuth [deg]'),ylabel('elevation [deg]'),zlabel('array factor [dB]')
    
    % 绘制水平/垂直方向切片
    subplot(3,2,[1 2])
    ind = find(elq==0);
    plot(azq(ind),20*log10(abs(afMat(ind)))), grid on
    set(gca,'XTick',[-135:45:135]), xlim([-135 135])
    xlabel('azimuth [deg]'),ylabel('array factor [dB]')
    title(sprintf('Rx 阵列波束，水平切片'))
    
    subplot(3,2,[3 4])
    ind = find(azq==0);
    plot(elq(ind),20*log10(abs(afMat(ind)))), grid on
    xlabel('elevation [deg]'),ylabel('array factor [dB]')
    title(sprintf('Rx 阵列波束，垂直切片'))   
    
    % 绘制三维阵列几何
    subplot(3,2,6)
    plot3(rx.acoord(:,1)*1E3,rx.acoord(:,2)*1E3,rx.acoord(:,3)*1E3,'o')
    grid on, axis square equal
    xlabel('X [mm]'),ylabel('Y [mm]'),zlabel('Z [mm]')
    title(sprintf('Rx 阵列几何, 共 %d 元素, 类型: %s',Nant,par.rx.antennaType))
end



%% 生成信道系数（CIRs）%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 移除“掉线/失效”的链路（Outage）：某些链路没有有效的路径数 Npaths（为 NaN）
OutageIdx = cellfun(@(x) isnan(x.Npaths), {pathData{par.linkSet}});
par.linkSet = par.linkSet(~OutageIdx);      % 仅保留非掉线链路
Nlink = length(par.linkSet);                % 实际参与仿真的链路数

% 统计每条链路的“tap数量”（即可解析的多径条数），以便统一预分配数组尺寸
for i_link = 1:Nlink
    tmpTaps(i_link) = length(pathData{par.linkSet(i_link)}.specPaths_UsedDyR.delay);
end
maxTap = max(tmpTaps);  % 所有链路中最大的 tap 数（统一维度）

% 预分配 CIR 与 Delay 的存储数组
% CIR 维度： [链路数 × Rx天线数 × Tx天线数 × 多径条数 × 时间样本数]
CIR   = zeros(Nlink, rx.Nant, tx.Nant, maxTap, par.Ntime);
% Delay 维度： [链路数 × 多径条数 × 时间样本数]（本模型中延迟不随时间变）
Delay = zeros(Nlink, maxTap,par.Ntime);

% 遍历每条链路，生成这条链路的多径参数与天线增益
for i_link = 1:Nlink
    % 取出该链路的路径功率/时延/角度等“测量得到的双向参数”
    powe = pathData{par.linkSet(i_link)}.specPaths_UsedDyR.power;    % 功率[dB]
    dela = pathData{par.linkSet(i_link)}.specPaths_UsedDyR.delay;    % 时延[ns]（注意单位）
    aoa  = pathData{par.linkSet(i_link)}.specPaths_UsedDyR.RxPhi;    % 到达方位角 AoA [deg]
    zoa  = pathData{par.linkSet(i_link)}.specPaths_UsedDyR.RxTheta;  % 到达“theta”角（0~360）
    aod  = pathData{par.linkSet(i_link)}.specPaths_UsedDyR.TxPhi;    % 发射方位角 AoD [deg]
    zod  = pathData{par.linkSet(i_link)}.specPaths_UsedDyR.TxTheta;  % 发射“theta”角（0~360）

    % 将原始的 theta(0~360°) 换算为传统“俯仰角”（Elevation，范围 -90°~+90°）
    % 约定：zoa=0 指向“上方”（以 90° 为水平基准进行对称映射）
    eoa = zoa;
    eoa(zoa>=0   & zoa<=180) = 90 - zoa(zoa>=0   & zoa<=180);
    eoa(zoa>180  & zoa<=360) =      zoa(zoa>180  & zoa<=360) - 270;
    eod = zod;
    eod(zod>=0   & zod<=180) = 90 - zod(zod>=0   & zod<=180);
    eod(zod>180  & zod<=360) =      zod(zod>180  & zod<=360) - 270;

    % 当前链路的有效多径条数（也可用 length(dela)）
    Npath = pathData{par.linkSet(i_link)}.Npaths;

    % 保存该链路各径的时延（原始单位为ns
    Delay(i_link, 1:Npath) = dela;

    %============================================= 
    % 初始光程d0与运动方向u的计算
    % 取初始 Tx/Rx 坐标
    txc  = pathData{par.linkSet(i_link)}.TxCoord(:).';   % [1x3]
    rxc0 = pathData{par.linkSet(i_link)}.RxCoord(:).';   % [1x3]
    
    % 运动方向的单位向量（与 par.DoTaz / par.DoTel 一致）
    ux = cosd(par.DoTaz)*cosd(par.DoTel);
    uy = sind(par.DoTaz)*cosd(par.DoTel);
    uz = sind(par.DoTel);
    u  = [ux, uy, uz];   % [1x3]
    
    % 每条径的初始光程 d0：如果有散射点，用 Tx->TxScat + RxScat->Rx，否则按 LOS 用 Tx->Rx
    d0 = zeros(1, Npath);
    TxSc = pathData{par.linkSet(i_link)}.specPaths_UsedDyR.TxScatLoc; % [3 x Npath] or []
    RxSc = pathData{par.linkSet(i_link)}.specPaths_UsedDyR.RxScatLoc; % [3 x Npath] or []
    
    for i_path = 1:Npath
        hasSc = ~isempty(TxSc) && ~isempty(RxSc) && ~isempty(TxSc(:,i_path));
        if hasSc
            txs = TxSc(:,i_path).';
            rxs = RxSc(:,i_path).';
            d0(i_path) = norm(txc - txs) + norm(rxs - rxc0);
        else
            d0(i_path) = norm(txc - rxc0);  % 视作 LOS
        end
    end
    %=======================================================

    % 为每条路径生成一个随机初相位（[0,1) 均匀分布），形成复指数相量
    phasor = exp(1j * 2*pi * rand(1, Npath));

    % 找到功率最大的路径索引，用于波束朝向“自动对准最大功率方向（orientA='max'）”
    [~, ind] = max(powe);
    % 记下 Tx/Rx 对应该最强径的（俯仰，方位）角度，四舍五入到整数度用于矩阵平移
    orientTx(i_link, :) = round([eod(ind) aod(ind)]);
    orientRx(i_link, :) = round([eoa(ind) aoa(ind)]);

    %% 计算发射端天线增益/导向矢量 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmpi(par.tx.antennaType, 'single3GPP')
        % 单天线（3GPP方向图）
        if strcmpi(tx.orientA, 'max')
            % 将方向图矩阵平移，使其主瓣指向最强路径方向（二维 circshift）
            A_dB_rot = circshift(tx.A_dB, orientTx(i_link, :));
        end
        % 插值：在（AoD, EoD）处采样旋转后的方向图，得到每条路径的天线增益[dB]
        antGain_dB  = interp2(az, el, A_dB_rot, aod, eod);
        % 转为“幅度增益”（线性幅度）
        tx.antGain  = 10.^(antGain_dB/20);
        % 记录当前链路下（“旋转后”的）天线坐标（单天线此处仅为占位）
        tx.rotAcoord(:, :, i_link) = tx.acoord;

    else
        % 阵列（URA/UCA）：使用几何坐标生成每条路径对应的导向矢量/阵列响应
        if strcmpi(tx.orientA, 'max')
            % 将阵列整体旋转，使其“波束朝向”最强路径方向（先绕 y，再绕 z）
            rotEle = roty(orientTx(i_link, 1));
            rotAzi = rotz(orientTx(i_link, 2));
            tmp    = rotEle * tx.acoord.';   % elevation 旋转
            tmp    = rotAzi * tmp;           % azimuth   旋转
            acoord = tmp.';                  % 旋转后的阵元坐标
        else
            % 保持用户指定的固定朝向（或 UCA 不用旋转）
            acoord = tx.acoord;
        end

        % 针对每一条路径，计算该方向上的导向矢量（阵列响应）
        for i_path = 1:Npath
            angAzSt = aod(i_path);  % 该径的发射方位角
            angElSt = eod(i_path);  % 该径的发射俯仰角
            % 该方向对应的“波矢量”（单位：rad/m），见平面波模型 k·r 相位
            gamma = -[ cosd(angAzSt).*cosd(angElSt), ...
                       sind(angAzSt).*cosd(angElSt), ...
                       sind(angElSt) ] * 2*pi/lambda;
            % 阵列响应：对每个阵元坐标 r，计算 e^{j * k·r}
            vecSt = exp(1j * gamma * acoord.');
            % 保存：Tx 端阵列增益/导向矢量（大小：Nant × Npath，复数）
            tx.antGain(:, i_path) = vecSt.';
        end
        % 记录本链路下旋转后的Tx阵列坐标
        tx.rotAcoord(:, :, i_link) = acoord;
    end

    %% 计算接收端天线增益/导向矢量 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmpi(par.rx.antennaType, 'single3GPP')
        % 单天线（3GPP方向图）
        if strcmpi(rx.orientA, 'max')
            % 同理：朝向最强路径
            A_dB_rot = circshift(rx.A_dB, orientRx(i_link, :));
        end
        % 在（AoA, EoA）处采样旋转后的方向图
        antGain_dB = interp2(az, el, A_dB_rot, aoa, eoa);
        rx.antGain = 10.^(antGain_dB/20);   % 幅度增益（线性）
        rx.rotAcoord(:, :, i_link) = rx.acoord;

    else
        % 阵列（URA/UCA）：同 Tx，按路径方向生成接收端导向矢量
        if strcmpi(rx.orientA, 'max')
            rotEle = roty(orientRx(i_link, 1));   % elevation 旋转
            rotAzi = rotz(orientRx(i_link, 2));   % azimuth   旋转
            tmp    = rotEle * rx.acoord.';
            tmp    = rotAzi * tmp;
            acoord = tmp.';
        else
            acoord = rx.acoord;
        end

        for i_path = 1:Npath
            angAzSt = aoa(i_path);  % 该径的到达方位角
            angElSt = eoa(i_path);  % 该径的到达俯仰角
            gamma = -[ cosd(angAzSt).*cosd(angElSt), ...
                       sind(angAzSt).*cosd(angElSt), ...
                       sind(angElSt) ] * 2*pi/lambda;
            vecSt = exp(1j * gamma * acoord.');   % 接收端导向矢量
            rx.antGain(:, i_path) = vecSt.';      % 保存：Nant × Npath
        end
        rx.rotAcoord(:, :, i_link) = acoord;      % 记录旋转后坐标
    end

    %% 根据“模型类型”生成 CIR（是否包含多普勒/时间演化）%%%%%%%%%%%%%%%%%%%%
    switch par.modelType
        case 'snapshot'  % 静态快照：每条链路仅 1 个时刻样本
            % CIR 维度： [link, rxAnt, txAnt, path, time]
            % 公式：h = a_rx(θ) * ( |α| * e^{jφ} ) * a_tx(θ)^T
            for i_path = 1:Npath
                CIR(i_link, :, :, i_path, 1) = ...
                    rx.antGain(:, i_path) * ...                      % Rx 导向矢量（列）
                    (db2mag(powe(i_path)) .* phasor(i_path)) * ...   % 路径复增益（幅度由测量功率，初相随机）
                    tx.antGain(:, i_path).';                         % Tx 导向矢量（行）
            end

        case 'virtualMotion'  % 虚拟运动：固定空间参数 + 多普勒/时间起伏
            % 每条路径的多普勒：fd = v/λ * cos(到达方向与运动方向的夹角)
            % 这里拆成方位/俯仰两个夹角余弦：cos(DoTaz-aoa)*cos(DoTel-eoa)
            fdop = cosd(par.DoTaz - aoa) .* cosd(par.DoTel - eoa) ...
                   * par.speed / SoL * par.fc;    % [Hz]

            % 采样间隔：2*SD 倍过采样（Nyquist 上的安全裕度）
            DeltaT = 1 / (2 * SD * par.speed / SoL * par.fc);

            % ========================================
            % 修改，增大可见衰落幅度

            % 选主径（ind = argmax power）
            d0_main = d0(ind);  % 循环前已算好的每径初始光程 d0
            cospsi  = abs(cosd(par.DoTaz - aoa(ind)) * cosd(par.DoTel - eoa(ind))); % 与主径夹角
            
            % 设定要达到的起伏幅度（dB）
            delta_dB = 1.0;  
            
            % 计算需要的样本数 Ntime
            Ntime_req = ceil( (2*SD*d0_main/(lambda*max(cospsi,1e-6))) * (10^(delta_dB/20) - 1) );
        
            if Ntime_req > par.Ntime
                par.Ntime = Ntime_req;
            end

            tvec = (0:par.Ntime-1) * DeltaT;  % 用更新后的 Ntime 重建时间向量
            % tvec   = (0:par.Ntime-1) * DeltaT;    % 时间轴（秒）

            % % 公式：在 snapshot 基础上乘以多普勒相位项 e^{j2π f_d t}
            % for i_path = 1:Npath
            %     for i_time = 1:par.Ntime
            %         CIR(i_link, :, :, i_path, i_time) = ...
            %             rx.antGain(:, i_path) * ...
            %             (db2mag(powe(i_path)) .* phasor(i_path)) * ...
            %             tx.antGain(:, i_path).' * ...
            %             exp(1j * 2*pi * fdop(i_path) * tvec(i_time));
            %     end
            % end
                
            % 修改=== 几何光程驱动的 FSPL 幅度 + 几何相位（隐含多普勒） ===
            for i_time = 1:par.Ntime
                tnow = tvec(i_time);
        
                % RX 位置随时间更新（目前TX不动）
                rxc_t = rxc0 + par.speed * tnow * u;   % [1x3]
        
                for i_path = 1:Npath
                    % 计算当前光程 d_t
                    hasSc = ~isempty(TxSc) && ~isempty(RxSc) && size(TxSc,2)>=i_path && ~isempty(TxSc(:,i_path));
                    if hasSc
                        txs = TxSc(:,i_path).';
                        rxs = RxSc(:,i_path).';
                        d_t = norm(txc - txs) + norm(rxs - rxc_t);   % 单反射或测得散射点
                    else
                        d_t = norm(txc - rxc_t);                     % LOS
                    end
        
                    % % 幅度变化（相对 FSPL）：|α(t)| = |α_meas| * (d0/d_t)
                    % amp_scale = d0(i_path) / max(d_t, 1e-6);  % 防止除零  幅度相位因子=d0/dt
                    % amp_t     = db2mag(powe(i_path)) * amp_scale; % 单径时变幅度 =根号p·幅度相位因子
                    % % 相位变化：phi(t) = phi0 - 2π fc (d_t - d0) / c
                    % phase_t = phasor(i_path) * exp(-1j * 2*pi * par.fc * (d_t - d0(i_path)) / SoL);
                    % % 叠加单径贡献（不更新导向矢量）
                    % CIR(i_link, :, :, i_path, i_time) = ...
                    %     rx.antGain(:, i_path) * (amp_t * phase_t) * tx.antGain(:, i_path).';


                    % 在距离过近时CIR数据产生了NaN值，因此设置与散射点的距离下限，避免奇点
                    d0(i_path) = max(d0(i_path), par.d_min);
                    if hasSc
                        d2_0 = max(norm(rxs - rxc0), par.d_min);
                        d2_t = max(norm(rxs - rxc_t), par.d_min);
                        d1 = norm(txc - txs);
                        d_t  = d1 + d2_t;  
                        amp_scale = d2_0 / d2_t;
                    else
                        d_t  = max(norm(txc - rxc_t), par.d_min);
                        amp_scale = d0(i_path) / d_t;
                    end
                    % 幅度限幅
                    amp_scale = min(max(amp_scale, par.amp_floor), par.amp_clip);
                    
                    % 复增益
                    amp_t   = db2mag(powe(i_path)) * amp_scale;
                    phase_t = phasor(i_path) * exp(-1j * 2*pi * par.fc * (d_t - d0(i_path)) / SoL);
                    
                    g_rx = rx.antGain(:,i_path);
                    g_tx = tx.antGain(:,i_path);
                    
                    if ~all(isfinite([amp_t, phase_t])) || any(~isfinite(g_rx)) || any(~isfinite(g_tx))
                        CIR(i_link,:,:,i_path,i_time) = 0;     % 兜底
                    else
                        CIR(i_link,:,:,i_path,i_time) = g_rx * (amp_t * phase_t) * g_tx.';
                    end
                end
            end

        case 'continuous'
            % 预留：如果未来做“链路间插值”的动态信道，可在此实现
            % 当前未实现
    end % end of modelType
end % end of links loop


%% Transform to frequency domain 
% 将时域信道脉冲响应 (CIR) 转换到频域信道频率响应 (CFR)

Deltaf = 5E6;         % 频点间隔 [Hz]，这里设置为 0.5 MHz 5mhz
freq = [par.fc-par.BW/2:Deltaf:par.fc+par.BW/2]; 
                        % 频率采样点集合 [Hz]
                        % 从中心频率 par.fc - 带宽的一半到 par.fc + 带宽的一半
                        % 间隔为 Deltaf
Nfreq = length(freq);   % 频率点的数量

% 初始化 CFR 矩阵
% 维度：链路数 × 接收天线数 × 发射天线数 × 频点数 × 时间点数
CFR = zeros(Nlink,rx.Nant,tx.Nant,Nfreq,par.Ntime);

% 遍历每条链路
for i_link = 1:Nlink
    tau = Delay(i_link,:,1)';     % 当前链路的路径时延向量 [s] (假设时延不随时间变化)

    % 构造指数矩阵，用于做 DFT 转换
    % exp(-j*2*pi*f*tau)，相当于将 CIR 投影到频域
    expmat = exp(-j*2*pi*tau*freq);   

    fprintf('Performing DFT %d/%d\n',i_link,Nlink)  % 输出进度信息

    % 遍历时间点
    for i_t = 1:par.Ntime
        % 遍历接收天线
        for i_rx = 1:rx.Nant
            % 遍历发射天线
            for i_tx = 1:tx.Nant
                % 提取 CIR(i_link, i_rx, i_tx, :, i_t)，维度是 Npath
                % 乘以 expmat 得到频域响应
                % 结果是 Nfreq 个频率点的 CFR
                CFR(i_link,i_rx,i_tx,:,i_t) = ...
                    squeeze(CIR(i_link,i_rx,i_tx,:,i_t)).' * expmat;  
            end % FOR tx
        end % FOR rx
    end % FOR time
end

% ============================================================
% 注：如果路径时延随时间变化（continuous 模式），
% 就需要在每个时间点重新生成 exp(-j*2*pi*f*tau)，
% 否则默认假设 Delay 不变，只计算一次。
% ============================================================

% % % 如果路径时延随时间变化，则使用下面这个较慢的循环
% % if strcmpi(par.modelType,'continuous')
% %     for i_link = 1:Nlink
% %         for i_t = 1:par.Ntime
% %             tau = Delay(i_link,:,i_t)'*1E9;     % 时延向量 [s]
% %             expmat = exp(-j*2*pi*tau*freq);     % 指数矩阵
% %             CFR(i_link,:,:,i_t) = squeeze(CIR(i_link,:,:,i_t))*expmat;  
% %             % 得到该时间点的频域信道响应
% %         end
% %     end
% % end


%% Visualization   
% 可视化信道模型结果
if visualizeChModel
    linkInd = 1;  % 选择要可视化的链路（对应 par.linkSet 中的索引）

    % ================= 功率延迟谱（PDP） =================
    figure(5),clf,hold on
    for i_link = 1:Nlink
        % 绘制每条链路的路径延迟与增益
        stem(Delay(i_link,:,1)*1E9, ... % 横轴：延迟 [ns]
             squeeze(20*log10(abs(CIR(i_link,1,1,:,1))))', ... % 纵轴：增益[dB]
             'BaseValue',-180) % 基线设为 -180 dB
        leg{i_link} = sprintf('Link %d',par.linkSet(i_link)); % 图例
    end
    grid on, hold off, xlim([0 max(get(gca,'XLim'))])
    xlabel('path delay [ns]'),ylabel('path gain [dB]'), 
    title('PDP, single antenna')
    legend(leg)
    
    % ================= 空间域可视化（角度-增益分布） =================
    figure(6),clf,hold on
    for i_link = 1:Nlink
        % 提取路径角度参数
        aoa  = pathData{par.linkSet(i_link)}.specPaths_UsedDyR.RxPhi;    % 到达方位角 [deg]
        zoa  = pathData{par.linkSet(i_link)}.specPaths_UsedDyR.RxTheta;  % 到达仰角 [deg]
        aod  = pathData{par.linkSet(i_link)}.specPaths_UsedDyR.TxPhi;    % 发射方位角 [deg]
        zod  = pathData{par.linkSet(i_link)}.specPaths_UsedDyR.TxTheta;  % 发射仰角 [deg]

        % 将 [0,360] 的仰角转换为 [-90,90] 的俯仰角
        eoa = zoa;  % elevation of arrival
        eoa(zoa>=0 & zoa<=180) = 90-zoa(zoa>=0 & zoa<=180);
        eoa(zoa>180 & zoa<=360) = zoa(zoa>180 & zoa<=360)-270;
        eod = zod;  % elevation of departure
        eod(zod>=0 & zod<=180) = 90-zod(zod>=0 & zod<=180);
        eod(zod>180 & zod<=360) = zod(zod>180 & zod<=360)-270;

        % 当前链路的路径数
        Npath = pathData{par.linkSet(i_link)}.Npaths; 
        leg{i_link} = sprintf('Link %d',par.linkSet(i_link));

        % 绘制 3D 棒状图，展示不同角度组合下的增益
        subplot(2,2,1),hold on
        stem3(aoa,aod,squeeze(20*log10(abs(CIR(i_link,1,1,1:Npath,1))))','BaseValue',-180)
        xlabel('AoA [deg]'),ylabel('AoD [deg]'),zlabel('gain [dB]'), grid on, view([20 60])

        subplot(2,2,2),hold on
        stem3(eoa,eod,squeeze(20*log10(abs(CIR(i_link,1,1,1:Npath,1))))','BaseValue',-180)
        xlabel('EoA [deg]'),ylabel('EoD [deg]'),zlabel('gain [dB]'), grid on, view([20 60])

        subplot(2,2,3),hold on
        stem3(aoa,eoa,squeeze(20*log10(abs(CIR(i_link,1,1,1:Npath,1))))','BaseValue',-180)
        xlabel('AoA [deg]'),ylabel('EoA [deg]'),zlabel('gain [dB]'), grid on, view([20 60])

        subplot(2,2,4),hold on
        stem3(aod,eod,squeeze(20*log10(abs(CIR(i_link,1,1,1:Npath,1))))','BaseValue',-180)
        xlabel('AoD [deg]'),ylabel('EoD [deg]'),zlabel('gain [dB]'), grid on
    end
    legend(leg), view([20 60])
    
    % ================= 如果是虚拟运动模型，显示时变结果 =================
    if strcmpi(par.modelType,'virtualMotion') 
        % CIR 随时间的变化
        figure(7)
% =================修改：使横轴为时间轴s而不是样本数
        plot(20*log10(abs(squeeze(CIR(linkInd,1,1,:,:))))')
        xlabel('time [s]'),ylabel('path gain [dB]'), 
        plot(tvec, 20*log10(abs(squeeze(CIR(linkInd,1,1,:,:))).'))
        xlabel('time [s]'),ylabel('path gain [dB]'), 
        title(sprintf('Link %d, single antenna, all paths',par.linkSet(linkInd)))
        
        % 所有链路在中心频率处的增益随时间变化
        figure(8)
        plot(20*log10(abs(squeeze(sum(CIR(:,1,1,:,:),4))))')
        xlabel('time [s]'),ylabel('gain at f_c [dB]'), title('All links, single antenna, @f_c')
        legend(leg)

        % CFR 随时间和频率的二维分布
        figure(9)
        surf(freq*1E-9,tvec,20*log10(abs(squeeze(CFR(1,1,1,:,:))))','EdgeColor','flat')
        xlabel('frequency [GHz]'),ylabel('time [s]'),zlabel('magnitude [dB]')
        title(sprintf('Link %d, single antenna',par.linkSet(linkInd)))

        % Doppler-频率二维谱
        figure(10)
        surf(linspace(-SD/2,SD/2,par.Ntime),freq*1E-9,...
            20*log10(abs(fftshift(fft(squeeze(CFR(linkInd,1,1,:,:)),[],2)))),'EdgeColor','flat')
        xlabel('normalized Doppler freq.'),ylabel('frequency [GHz]')
        zlabel('magnitude [dB]')
        title(sprintf('Link %d, single antenna',par.linkSet(linkInd)))
    end

    % ================= 如果是快照模型，显示 CFR 结果 =================
    if strcmpi(par.modelType,'snapshot') 
        figure(14)
        plot(freq*1E-9,20*log10(abs(squeeze(CFR(:,1,1,:,1)))))
        xlabel('frequency [GHz]'),ylabel('magnitude [dB]')
        title('CFR, single time instant, single antenna')
        grid on, legend(leg)
    end
end 

%% Visualize layout
% 可视化 Tx/Rx/散射体的空间布局
if visualizeLayout
    for i_link = 1:Nlink
        figNo = 100 + i_link;
        Npath = pathData{par.linkSet(i_link)}.Npaths;
        txc  = pathData{par.linkSet(i_link)}.TxCoord;   % Tx 坐标 [m]
        rxc  = pathData{par.linkSet(i_link)}.RxCoord;   % Rx 坐标 [m]

        % ========== 绘制 Tx、Rx 和散射体位置 ==========
        figure(figNo)
        plot3(txc(1),txc(2),txc(3),'o','MarkerFaceColor','r','MarkerSize',10), hold on
        plot3(rxc(1),rxc(2),rxc(3),'^','MarkerFaceColor','g','MarkerSize',10)

        % 绘制 Tx/Rx 方向向量
        distTR = sqrt(sum((txc-rxc).^2));   % Tx-Rx 距离
        [otx,oty,otz] = sph2cart(deg2rad(orientTx(i_link,2)), deg2rad(orientTx(i_link,1)), distTR/8);
        [orx,ory,orz] = sph2cart(deg2rad(orientRx(i_link,2)), deg2rad(orientRx(i_link,1)), distTR/8);
        plot3([txc(1) txc(1)+otx], [txc(2) txc(2)+oty], [txc(3) txc(3)+otz], '-','LineWidth',2)
        plot3([rxc(1) rxc(1)+orx], [rxc(2) rxc(2)+ory], [rxc(3) rxc(3)+orz], '-','LineWidth',2)

        % 绘制散射点位置
        Nscat = 0;
        for i_path = 1:Npath
            txSc = pathData{par.linkSet(i_link)}.specPaths_UsedDyR.TxScatLoc(:,i_path);
            rxSc = pathData{par.linkSet(i_link)}.specPaths_UsedDyR.RxScatLoc(:,i_path);
            if ~isempty(txSc)
                Nscat = Nscat+1;    % 有效散射点计数
                plot3(txSc(1),txSc(2),txSc(3),'o', rxSc(1),rxSc(2),rxSc(3),'^');
                plot3([txc(1), txSc(1)],[txc(2), txSc(2)],[txc(3), txSc(3)],'k:',...
                      [rxc(1), rxSc(1)],[rxc(2), rxSc(2)],[rxc(3), rxSc(3)],'k--')
            end
        end
        hold off, grid on, axis equal
        title(sprintf('Link %d, scatterers located %d/%d',par.linkSet(i_link),Nscat,Npath))
        legend('Tx','Rx','Tx orientation','Rx orientation','Tx scat','Rx scat', 'Location','Best')
        xlabel('X [m]'),ylabel('Y [m]'),zlabel('Z [m]')

        % ========== 绘制 Tx/Rx 阵列布局 ==========
        figure(figNo+100),
        % Tx 子图
        subplot(2,1,1)
        if strcmp(par.tx.antennaType,'single3GPP')
            plot3(txc(1),txc(2),txc(3),'o','MarkerFaceColor','r','MarkerSize',10),hold on
        else % URA/UCA
            plot3(txc(1)+tx.rotAcoord(:,1,i_link),txc(2)+tx.rotAcoord(:,2,i_link),...
                  txc(3)+tx.rotAcoord(:,3,i_link),'o','MarkerSize',4),hold on
        end
        plot3([txc(1) txc(1)+otx*0.1], [txc(2) txc(2)+oty*0.1], [txc(3) txc(3)+otz*0.1], '-','LineWidth',2)
        hold off, axis equal, grid on
        xlabel('X [m]'),ylabel('Y [m]'),zlabel('Z [m]')
        title(sprintf('Link %d, Tx: %s',par.linkSet(i_link),par.tx.antennaType))
        legend('Tx','Tx orientation','Location','Best')

        % Rx 子图
        subplot(2,1,2)
        if strcmp(par.rx.antennaType,'single3GPP')
            plot3(rxc(1),rxc(2),rxc(3),'^','MarkerFaceColor','g','MarkerSize',10),hold on
        else % URA/UCA
            plot3(rxc(1)+rx.rotAcoord(:,1,i_link),rxc(2)+rx.rotAcoord(:,2,i_link),...
                  rxc(3)+rx.rotAcoord(:,3,i_link),'^','MarkerSize',4),hold on
        end
        plot3([rxc(1) rxc(1)+orx*0.1], [rxc(2) rxc(2)+ory*0.1], [rxc(3) rxc(3)+orz*0.1], '-','LineWidth',2)
        hold off, axis equal, grid on
        xlabel('X [m]'),ylabel('Y [m]'),zlabel('Z [m]')
        title(sprintf('Link %d, Rx: %s',par.linkSet(i_link),par.rx.antennaType))
        legend('Rx','Rx orientation', 'Location','Best')
    end
end


