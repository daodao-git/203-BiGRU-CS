function [RxScatLoc,TxScatLoc] = randLoc(RxCoord,TxCoord,RxPhi,TxPhi,RxTheta,TxTheta,sRx,sTx,del)
%s - scale parameter of rayleigh distribution
%phi - azimuth angle in degrees
%theta - polar angle in degrees
%Assumes double bounce scattering

pd = makedist('Rayleigh','B',sRx);
d = random(pd,1,length(RxPhi));
[rx,ry,rz] = sph2cart(RxPhi*(pi/180),(90-RxTheta)*(pi/180),d);
r = [rx;ry;rz];
RxScatLoc = RxCoord+r; %Initial random RxScatLoc

pd = makedist('Rayleigh','B',sTx);
d = random(pd,1,length(TxPhi));
[rx,ry,rz] = sph2cart(TxPhi*(pi/180),(90-TxTheta)*(pi/180),d);
r = [rx;ry;rz];
TxScatLoc = TxCoord+r; %Initial random TxScatLoc

%Adjust Tx and Rx ScatLocs so that the overall path delay is equal to the
%measured delay del.
RxScatVec = RxScatLoc - RxCoord;
TxScatVec = TxScatLoc - TxCoord;
RxScatu = RxScatVec./norm(RxScatVec);
TxScatu = TxScatVec./norm(TxScatVec);
%TxScatLoc is adjusted while RxScatLoc is held stationary
testTxScat = TxCoord + TxScatu*[0.5:0.1:100]; %equally spaced points along the line formed between TxCoord and 100meters away from Tx along TxScatu 
XcsDist = vecnorm(RxScatVec) + vecnorm(testTxScat - TxCoord) + vecnorm(testTxScat-RxScatLoc) - 299792458*del;
minXcsDist = min(abs(XcsDist));
ind = find(abs(XcsDist) == minXcsDist);
TxScatLoc = testTxScat(:,ind); % Adjusted TxScatLoc
if minXcsDist > 0.1 %If the minimum excess distance by moving TxScatLoc is above some threshold, RxScatLoc is then adjusted next.
    TxScatVec = TxScatLoc - TxCoord;
    testRxScat = RxCoord + RxScatu*[0.5:0.1:100];
    XcsDist2 = vecnorm(testRxScat - RxCoord) + vecnorm(TxScatVec) + vecnorm(testRxScat-TxScatLoc) - 299792458*del;   
    RxScatLoc = testRxScat(:,find(abs(XcsDist2) == min(abs(XcsDist2))));
end
end