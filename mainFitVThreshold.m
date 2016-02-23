function [] = mainFitVThreshold(cName)

% 1. INITIALIZATION
dataName = sprintf('GIFTraining_%s',cName);
load(dataName)

tic

samplingFreq = GIF.samplingFreq;
dt = 1e3/samplingFreq;
tRefr = round(GIF.tRefr/dt);

TGamma = GIF.TGamma; nbrBinGamma = GIF.nbrBinGamma;

spike = GIF.spike;
[S,GammaPattern,tGammaBin] = buildSMatrix(nbrBinGamma,TGamma,spike,tRefr*dt,'gammaFS',samplingFreq);
tGamma = 0:dt:length(GammaPattern)*dt - dt;

firing_rate = full(sum(spike)/(length(GIF.v)/samplingFreq));
sprintf('firing rate %f [Hz]',firing_rate)

clear I

V = GIF.v;
figure(5),hold on,
t=1:1:length(V);
plot(t,V,'k'),plot(t(spike==1),V(spike==1),'.r')
spiketime = (find(spike==1));

for i=1:sum(spike)
    V(spiketime(i)+1:spiketime(i)+tRefr) = nan;
end
plot(t,V,'m')
ind = isnan(V); 
V(ind) = []; SS = S; SS(:,ind) = [];
spikeM = spike'; spikeM(ind) = [];

% 1. INITIALIZE WITH CST THRESHOLD

maxIter = 500;
L = nan(maxIter,1);
X = [V;ones(size(V))];
paramHat = [.15;11];
DV = 1/paramHat(1); V0 = -paramHat(2)*DV;
L(1) = likelihood(V',V0',spikeM,DV,dt);
lambda = .8; TolL = 1e-12;
iter=2;
while(iter<=maxIter)
    [gradL HL] = dLdtheta(X,paramHat,spikeM,dt);
    paramHat = paramHat - lambda*(HL\gradL);
    clear gradL HL
    
    DV = 1/paramHat(1);
    V0 = -paramHat(2)*DV;
    L(iter) = likelihood(V',V0',spikeM,DV,dt);
    DL = abs((L(iter-1)-L(iter))/L(iter-1));
    %     sprintf('L = %f, V0 = %f, DV = %f => DL = %f',L(iter),V0,DV,DL)
    
    if(DL<TolL)
        iter=maxIter+1;
    else
        iter = iter+1;
    end
end

paramHatCstThreshold = paramHat;
paramHat0 = [paramHatCstThreshold;zeros(nbrBinGamma,1)];

% 2. MOVING THRESHOLD
maxIter = 500; lambda = .8; TolL = 1e-12;

L = nan(maxIter,1);
X = [V;ones(size(V));SS];
paramHat = paramHat0;
DV = 1/paramHat(1); V0 = -paramHat(2)*DV;
gammaBin = -paramHat(3:end)'*DV;
L(1) = likelihood(V',(V0+(gammaBin*SS))',spikeM,DV,dt);
iter=2;
while(iter<=maxIter)
    [gradL,HL] = dLdtheta(X,paramHat,spikeM,dt);
    paramHat = paramHat - lambda*(HL\gradL);
    clear gradL HL
    
    DV = 1/paramHat(1);
    V0 = -paramHat(2)*DV;
    gammaBin = -paramHat(3:end)'*DV;
    L(iter) = likelihood(V',(V0+(gammaBin*SS))',spikeM,DV,dt);
    DL = abs((L(iter-1)-L(iter))/L(iter-1));
    
    sprintf('L = %f, V0 = %f, DV = %f => DL = %f',L(iter),V0,DV,DL)
    figure(6),hold on
    subplot(2,2,1),hold on,plot(iter,L(iter),'.k'),ylabel('L')
    subplot(2,2,2),hold on,plot(iter,V0,'.k'),ylabel('V0')
    subplot(2,2,3),hold on,plot(iter,DV,'.k'),ylabel('DV')
    subplot(2,2,4),hold on,plot(tGamma,gammaBin*GammaPattern,'k'),ylabel('Gamma')
    
    if(DL<TolL)
        iter=maxIter+1;
    else
        iter = iter+1;
    end
end

GIF.threshold.V0 = V0;
GIF.threshold.DV = DV; GIF.threshold.gammaExtracted = gammaBin*GammaPattern;
GIF.threshold.gammaBin = gammaBin; GIF.threshold.GammaPattern = GammaPattern;
GIF.threshold.tGamma = tGamma; GIF.threshold.tGammaBin = tGammaBin; GIF.threshold.L = L;
GIF.threshold.TGamma = TGamma; GIF.threshold.nbrBinGamma = nbrBinGamma;
GIF.threshold.gamma = GIF.threshold.gammaExtracted;
GIF.threshold.param = GIF.paramLIFEta;
GIF.threshold.param(6) = GIF.threshold.V0; GIF.threshold.param(7) = GIF.threshold.DV;
GIF.threshold.CPUTime = toc;

errorParam = 100*mean([abs((GIF.threshold.param - GIFRef.param)./GIFRef.param) ...
    abs((GIF.etaExtracted - GIFRef.eta)./GIFRef.eta) ...
    abs((GIF.threshold.gammaExtracted - GIFRef.gamma)./GIFRef.gamma)]);
GIF.eParam = errorParam;

fileName = sprintf('GIFTraining_%s',cName);
save(fileName,'GIF','GIFRef');

figure(7),hold on,
plot(GIFRef.tGamma,GIFRef.gamma,'k')
plot(GIF.threshold.tGammaBin,GIF.threshold.gammaBin,'.b'),plot(GIF.threshold.tGamma,GIF.threshold.gammaExtracted,'b')
ylabel('gamma [mV]'),xlabel('time [ms]')

end