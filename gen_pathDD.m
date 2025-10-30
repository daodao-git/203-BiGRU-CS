function [pathsData] = gen_pathDD(scenario)
%OUTPUT DATA "pathsData" - see the last part of this code for formats
addpath('functions\');

%LOADING DOUBLE-DIRECTIONAL PATH DATA FOR 140 GHZ ENTRANCE HALL
unpackStruct(GetPrm(scenario))
if ~valid
    disp('Scenario is not valid'); return;
end
disp(InitText)
load(dataDir)

pathsData = [];
for il = 1:length(RTdata) %link loop
    RxPhi = []; RxTheta = []; TxTheta = []; TxPhi = []; RxScatLoc = []; TxScatLoc = []; power = []; delay = []; NBounce = [];
    TxCoord = RTdata(il).TxCoord;
    RxCoord = RTdata(il).RxCoord;
    if RTdata(il).Npaths > 0
        bLOS = RTdata(il).bLOS;
        for ip = 1:RTdata(il).Npaths %path loop to generate TxPhi, TxTheta and RxTheta for un-indentified paths
            if RTdata(il).RTpaths(ip).NBounce == -1 %for paths not mapped on geometry
                %Rx direction seen from Tx
                Rxpos = RxCoord - TxCoord;
                [TxLOSaz,TxLOSel,~] = cart2sph(Rxpos(1),Rxpos(2),Rxpos(3)); %output is rad
                %Tx direction seen from Rx
                Txpos = TxCoord - RxCoord;
                [RxLOSaz,RxLOSel,~] = cart2sph(Txpos(1),Txpos(2),Txpos(3)); %output is rad
                
                del = RTdata(il).RTpaths(ip).del;

                %generate TxPhi, TxTheta and RxTheta
				RxPhiNew = RTdata(il).RTpaths(ip).ang;
				lap_rand = randl(1,3); %Laplacian random variable with zero-mean and one-std
                if bLOS
				    RxThetaNew = lap_rand(1)*bRxThetaLOS + muRxThetaLOS - rad2deg(RxLOSel);
				    TxPhiNew = lap_rand(2)*bTxPhiLOS + muTxPhiLOS + rad2deg(TxLOSaz);
				    TxThetaNew = lap_rand(3)*bTxThetaLOS + muTxThetaLOS - rad2deg(TxLOSel);
                    [RxScatLocNew,TxScatLocNew] = randLoc(RxCoord',TxCoord',RxPhiNew,TxPhiNew,RxThetaNew,TxThetaNew,sRxDistLOS,sTxDistLOS,del);
                else
				    RxThetaNew = lap_rand(1)*bRxThetaNLOS + muRxThetaNLOS - rad2deg(RxLOSel);
				    TxPhiNew = lap_rand(2)*bTxPhiNLOS + muTxPhiNLOS + rad2deg(TxLOSaz);
				    TxThetaNew = lap_rand(3)*bTxThetaNLOS + muTxThetaNLOS - rad2deg(TxLOSel);
                    [RxScatLocNew,TxScatLocNew] = randLoc(RxCoord',TxCoord',RxPhiNew,TxPhiNew,RxThetaNew,TxThetaNew,sRxDistNLOS,sTxDistNLOS,del);
                end

                RxPhi = [RxPhi RxPhiNew]; %Rx phi from measurements, deg
                RxTheta = [RxTheta RxThetaNew]; %Rx theta, deg
                TxPhi   = [TxPhi TxPhiNew]; %Tx phi, deg
                TxTheta = [TxTheta TxThetaNew]; %Tx theta, deg

                %generate RxScatLoc and TxScatLoc
                RxScatLoc = [RxScatLoc RxScatLocNew];
                TxScatLoc = [TxScatLoc TxScatLocNew];
            else
                RxPhi =   [RxPhi   RTdata(il).RTpaths(ip).RxPhi]; %Rx phi, deg
                RxTheta = [RxTheta RTdata(il).RTpaths(ip).RxTheta]; %Rx theta, deg
                TxPhi   = [TxPhi   RTdata(il).RTpaths(ip).TxPhi]; %Tx phi, deg
                TxTheta = [TxTheta RTdata(il).RTpaths(ip).TxTheta]; %Tx theta, deg
                RxScatLoc = [RxScatLoc RTdata(il).RTpaths(ip).Scat(:,1)];
                TxScatLoc = [TxScatLoc RTdata(il).RTpaths(ip).Scat(:,end)];
            end
        end
        power = [RTdata(il).RTpaths.gain]; %gain, dB
        delay = [RTdata(il).RTpaths.del]; %delay, ns
        NBounce = [RTdata(il).RTpaths.NBounce];
    end

    %populate data
    pathsData{il}.specPaths_UsedDyR.RxPhi   = RxPhi;                        % deg
    pathsData{il}.specPaths_UsedDyR.RxTheta = RxTheta;                      % polar angle (NOT elevation), deg
    pathsData{il}.specPaths_UsedDyR.TxPhi   = TxPhi;                        % deg
    pathsData{il}.specPaths_UsedDyR.TxTheta = TxTheta;                      % polar angle (NOT elevation), deg
    pathsData{il}.specPaths_UsedDyR.RxScatLoc = RxScatLoc;                  % last reflection point, m
    pathsData{il}.specPaths_UsedDyR.TxScatLoc = TxScatLoc;                  % first reflection point, m
    pathsData{il}.specPaths_UsedDyR.delay   = delay;                        % ns
    pathsData{il}.specPaths_UsedDyR.power   = power;                        % dB
    pathsData{il}.specPaths_UsedDyR.Nbounce = NBounce;
    pathsData{il}.scenario                  = RTdata(il).scenario;
    pathsData{il}.TxIdx                     = RTdata(il).TxIdx;
    pathsData{il}.TxName                    = RTdata(il).TxName;
    pathsData{il}.RxIdx                     = RTdata(il).RxIdx;
    pathsData{il}.bLOS                      = RTdata(il).bLOS;              % 1 for LOS, 0 for NLOS
    pathsData{il}.linkDist                  = RTdata(il).linkDist;          % meters
    pathsData{il}.TxCoord                   = RTdata(il).TxCoord;           % meters
    pathsData{il}.RxCoord                   = RTdata(il).RxCoord;           % meters
    pathsData{il}.DynamicRange              = RTdata(il).DynamicRange;      % dB
    pathsData{il}.OutageDynamicRange        = RTdata(il).OutageDynamicRange;
    pathsData{il}.UsedDynamicRange          = RTdata(il).UsedDynamicRange;
    pathsData{il}.fc                        = RTdata(il).fc;                % Hz
    pathsData{il}.BW                        = RTdata(il).BW;                % Hz
    pathsData{il}.Npaths                    = RTdata(il).Npaths;
end