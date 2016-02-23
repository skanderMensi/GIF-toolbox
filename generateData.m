function [] = generateData(cellName)
close('all'),clc

% PARAM REF FROM Pozzorini, Nature Neuroscience 2013
samplingFreq = 20000;
dt = 1e3/samplingFreq;

tau = 15.3;   %  ms
R = 93.2;    % MOhm
C = tau/R;  % nF
gL = 1/R;   % nS
EL = -69.4;   % mV
ER = -38.8;   % mV
tRefr = 4;  % ms

V0 = -51.9;   % mV
DV = 0.75;   % mV

% PARAM REF KERNEL
TEta = 5000;% 5 seconds
nbrBinEta = 26;% 26 bins

alphaEta = 0.44;   % nA
betaEta = 0.76;  %

spikeTemp = [1;zeros(TEta/dt,1)]';
[~,EtaPattern,tEtaBin] = buildSMatrix(nbrBinEta,TEta,spikeTemp,tRefr,'eta',samplingFreq); clear S

TGamma = 5000;% 5 seconds
nbrBinGamma = 26;% 26 bins
alphaGamma = 40;   % mV
betaGamma = 0.87;  %

spikeTemp = [1;zeros(TGamma/dt,1)]';
[~,GammaPattern,tGammaBin] = buildSMatrix(nbrBinGamma,TGamma,spikeTemp,tRefr,'gammaFS',samplingFreq); clear SS

% APPROXIMATE KERNEL WITH POWERLAW %
tEta = 0:dt:TEta-dt;
etaRef = alphaEta*tEta.^(-betaEta);
etaRef(1:20) = etaRef(20);
etaBinRefPL = etaRef(round(tEtaBin/dt));
etaRefPL = etaBinRefPL*EtaPattern;
loglog(tEta,etaRef,'k'),hold on,loglog(tEta,etaRefPL,'r')
loglog(tEtaBin,etaBinRefPL,'.r')

tGamma = 0:dt:TGamma-dt;
gammaRef = alphaGamma*tGamma.^(-betaGamma);
gammaRef(1:50) = gammaRef(20);
gammaBinRefPL = gammaRef(round(tGammaBin/dt));
gammaRefPL = gammaBinRefPL*GammaPattern;
loglog(tGamma,gammaRef,'k'),hold on,loglog(tGamma,gammaRefPL,'b')
loglog(tGammaBin,gammaBinRefPL,'.b')

GIFRef.C = C; GIFRef.gL = gL; GIFRef.tau = tau; GIFRef.R = R;
GIFRef.EL = EL; GIFRef.ER = ER; GIFRef.tRefr = tRefr;
GIFRef.V0 = V0; GIFRef.DV = DV;

GIFRef.eta = etaRefPL; GIFRef.tEta = tEta;
GIFRef.nbrBinEta = nbrBinEta; GIFRef.TEta = TEta;
GIFRef.gamma = gammaRefPL; GIFRef.tGamma = tGamma;
GIFRef.nbrBinGamma = nbrBinGamma; GIFRef.TGamma = TGamma;

GIFRef.param = [C gL EL ER tRefr V0 DV];

% GENERATE TRAINING SET @ ~10 Hz %

T = 100; % s

mu = 0.25;% nA
sigma = 0.25;% nA
dSigma = 0.5;
f = 0.2; % Hz
tau = 3; % ms

I = noisySinWave(mu,sigma,dSigma,f,tau,T,samplingFreq);
[spike,V] = IFEtaMTNu(GIFRef.param,GIFRef.eta',GIFRef.gamma',I,1,samplingFreq);
V = V';

rate = sum(spike)/(length(I)/samplingFreq);
tempDisp = sprintf('Training Set firing Rate: %.2f Hz',rate);
disp(tempDisp)

fileName = sprintf('%sTraining.mat',cellName);
save(fileName,'GIFRef','I','V','spike','samplingFreq');

% GENERATE TEST SET %

T = 10; % s
I = noisySinWave(mu,sigma,dSigma,f,tau,T,samplingFreq);
nbrRepet = 9;
V = nan(nbrRepet,T*samplingFreq);
spike = nan(nbrRepet,T*samplingFreq);
for i=1:nbrRepet
    [spike(i,:),V(i,:)] = IFEtaMTNu(GIFRef.param,GIFRef.eta',GIFRef.gamma',I,1,samplingFreq);
    rate = sum(spike(i,:))/(length(I)/samplingFreq);
    tempDisp = sprintf('Test Set %i firing Rate: %.2f Hz',i,rate);
    disp(tempDisp)
end

tTest = 0:1/samplingFreq:T - 1/samplingFreq;
figure(2),hold on,
subplot(2,1,1),hold on,plot(tTest,V(1,:),'k')
subplot(2,1,2),hold on,
for i=1:nbrRepet
    plot(tTest(spike(i,:)==1),i,'.k')
end

fileName = sprintf('%sTest.mat',cellName);
save(fileName,'GIFRef','I','V','spike','samplingFreq');

end