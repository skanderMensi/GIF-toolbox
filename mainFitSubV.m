function [] = mainFitSubV(cName)

% 0. INITIALIZATION
dataName = sprintf('%sTraining.mat',cName);
load(dataName)

tic

dt = 1e3/samplingFreq;
DeltaT = 4/dt;

tRefr = round(GIFRef.tRefr/dt);
VLimit = 0.0; % spike detection

TEta = GIFRef.TEta; nbrBinEta = GIFRef.nbrBinEta;

diffV = [0;diff(V')]/dt;
T = length(V)*dt;

[S,EtaPattern,tEtaBin] = buildSMatrix(nbrBinEta,TEta,spike,(tRefr*dt),'eta',samplingFreq);
tEta = 0:dt:length(EtaPattern)*dt - dt;

firing_rate = full(sum(spike)/(T*1e-3));
sprintf('firing rate %f [Hz]',firing_rate)

figure(3),hold on
t=1:1:length(V);
plot(t,V,'k')
plot(t(spike==1),V(spike==1),'.r')

% remove spikes
spiketime = (find(spike==1));
nbrSpike = length(spiketime);
Er = nan(nbrSpike,1);
for i=1:nbrSpike
    
    Er(i) = V(spiketime(i)+tRefr);
    V(spiketime(i)+1 - DeltaT:spiketime(i)+tRefr) = nan;
    
end
plot(t,V,'m')
ind = isnan(V);
V(ind) = []; diffV(ind) = []; I(ind) = [];
S(:,ind) = [];
Er = mean(Er);

sprintf('Er = %f, Initialization done in :',Er)

% 2. FIT LIF + SPIKE-TRIGGERED CURRENT ETA
X = [V; ones(size(V)); I; -1*S]';
[b,~,r] = regress(diffV,X); clear X bint
C = 1/b(3); gL = -b(1)*C; EL = (b(2)*C)/gL;
eta = (b(4:end)*C)'*EtaPattern;

sprintf('RMSE = %f => C = %f, gL = %f, EL = %f, tau = %f, in :',sqrt(mean(r.^2)),C,gL,EL,C/gL)

GIF.TEta = TEta; GIF.nbrBinEta = nbrBinEta;
GIF.TGamma = GIFRef.TGamma; GIF.nbrBinGamma = GIFRef.nbrBinGamma;
GIF.tRefr = tRefr*dt; GIF.C = C; GIF.gL = gL;
GIF.EL = EL; GIF.Er = Er; GIF.rmse = sqrt(mean(r.^2));
GIF.etaExtracted = eta; GIF.etaBin = (b(4:end)*C)';
GIF.EtaPattern = EtaPattern; GIF.tEta = tEta; GIF.tEtaBin = tEtaBin;
GIF.VLimit = VLimit;
GIF.eta = GIF.etaExtracted;

GIF.paramLIFEta = [C gL EL Er tRefr*dt]; GIF.etaLIFEta = GIF.eta';
GIF.cName = cName; GIF.samplingFreq = samplingFreq;

figure(4),hold on,
plot(GIFRef.tEta,GIFRef.eta,'k')
plot(GIF.tEtaBin,GIF.etaBin,'.r'),plot(GIF.tEta,GIF.etaExtracted,'r')
ylabel('eta [nA]'),xlabel('time [ms]')

GIF.CPUTime = toc;
GIF.CPUTime = GIF.CPUTime;
disp(GIF.CPUTime)

% 3 TEST MODELS ON TEST SET

tic
clear V I vIF vIFEta samplingFreq
dataName = sprintf('%sTest.mat',cName);
load(dataName)

VTest = V; ITest = I; clear I V 
nbrRepet = min(size(VTest)); T = max(size(VTest));
vIFEtaTest = nan(nbrRepet,T); rmseTestIFEta = nan(nbrRepet,1);
for i=1:nbrRepet
    
    spiketime = (find(spike(i,:)==1))';
    spiketime = spiketime - 1;
    
    vIFEtaTest(i,:) = IFEtaSpike(GIF.paramLIFEta,GIF.etaLIFEta,ITest(1,:),spiketime,samplingFreq);
    rmseTestIFEta(i) = sqrt(nanmean((VTest(i,:)-vIFEtaTest(i,:)).^2));
    
end
GIF.rmseTest = rmseTestIFEta;

clear VTest ITest vIFEtaTest rmseTestIFEta

dataName = sprintf('%sTraining.mat',cName);
load(dataName)

nbrRepet = min(size(I)); T = max(size(I));
vIFEta = nan(nbrRepet,T); rmseTrainingIFEta = nan(nbrRepet,1);

for i=1:nbrRepet

    spiketime = (find(spike(i,:)==1))';
    spiketime = spiketime - 1;
    
    vIFEta(i,:) = IFEtaSpike(GIF.paramLIFEta,GIF.etaLIFEta,I(i,:),spiketime,samplingFreq);
    rmseTrainingIFEta(i) = sqrt(nanmean((V(i,:)-vIFEta(i,:)).^2));
        
end
GIF.rmseTraining = rmseTrainingIFEta;

sprintf('training: rmseIFEta = %f mV, \n test: rmseIFEta = %f mV',...
    mean(GIF.rmseTraining),mean(GIF.rmseTest))

GIF.VData = V(1,:); GIF.IData = I(1,:); GIF.v = vIFEta(1,:); GIF.spike = spike(1,:);

fileName = sprintf('GIFTraining_%s',cName);
save(fileName,'GIF','GIFRef');

end