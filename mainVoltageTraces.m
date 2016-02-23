function [] = mainVoltageTraces(cName)
clc,close('all')

dataName = sprintf('%sTest.mat',cName);
load(dataName)
spikeData = spike;
fileName = sprintf('GIFTraining_%s',cName);
load(fileName)

nbrRepet = min(size(V));
spikeGIF=nan(size(V));
VGIF=nan(size(V));

for k=1:nbrRepet
    [spikeGIF(k,:),VGIF(k,:)] = IFEtaMTNu(GIF.threshold.param,GIF.eta',GIF.threshold.gamma',I(1,:),1,samplingFreq);
end

% Compute PSTH %
winSize = 0.5;
windowConv = ones(1,winSize*samplingFreq);
psthData = nan(min(size(spikeData)),length(windowConv)+length(spikeData)-1);
for i=1:min(size(spikeData))
    psthData(i,:) = conv(windowConv,spikeData(i,:))/winSize;
end

mPSTHData = mean(psthData,1); mPSTHData = mPSTHData(winSize*samplingFreq:end);
mPSTHGIF = conv(windowConv,mean(spikeGIF))/winSize;
mPSTHGIF = mPSTHGIF(winSize*samplingFreq:end);

t = 0:1/samplingFreq:length(V)/samplingFreq - 1/samplingFreq;
figure(1),hold on,
subplot(6,1,1),hold on,plot(t,I,'k'),ylabel('I (nA)')
subplot(6,1,2:3),hold on,plot(t,V,'k'),plot(t,VGIF(1,:),'r'),ylabel('Voltage (V)')
subplot(6,1,4:5),hold on
for i=1:nbrRepet
    plot(t(spikeData(i,:)==1),i,'.k')
    plot(t(spikeGIF(i,:)==1),-i,'.r')
end
subplot(6,1,6),hold on,plot(t,mPSTHData,'k'),plot(t,mPSTHGIF,'r')
ylabel('PSTH (Hz)'),xlabel('time (s)')

figure(2),hold on
subplot(1,2,1),hold on,plot(GIFRef.tEta,GIFRef.eta,'k')
plot(GIF.tEtaBin,GIF.etaBin,'.r'),plot(GIF.tEta,GIF.etaExtracted,'r')
ylabel('eta [nA]'),xlabel('time [ms]')
subplot(1,2,2),hold on,plot(GIFRef.tGamma,GIFRef.gamma,'k')
plot(GIF.threshold.tGammaBin,GIF.threshold.gammaBin,'.b'),plot(GIF.threshold.tGamma,GIF.threshold.gammaExtracted,'b')
ylabel('gamma [mV]'),xlabel('time [ms]')

end