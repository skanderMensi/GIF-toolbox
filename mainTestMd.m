function [] = mainTestMd(cName)

dataName = sprintf('%sTest.mat',cName);
load(dataName)

tic

dt = 1e3/samplingFreq;
delta = 4/dt;

% 1. COMPUTE MuDATA and NuDATA
mData = 0;
spikeData = spike;

nbrRepet = min(size(spikeData));
for i=1:nbrRepet
    for j=1:nbrRepet
        if(i<j)
            mData = mData+inprodGamma(full(spikeData(i,:)),full(spikeData(j,:)),delta);
        end
    end
end
mData = (2/((nbrRepet-1)*nbrRepet))*mData;
nuData = full(mean(spikeData,1));
sprintf('Data PSTH computed in: ')

% 2. GENERATE MODEL PSTH

fileName = sprintf('GIFTraining_%s',cName);
load(fileName);

% SETUP POOL
nbrRepet = 500;
% nbrRepet = 540;
% nbrPool = 12;
% nbrRepetPerPool = floor(nbrRepet/nbrPool);

tic
tempParam = GIF.threshold.param;
tempEta = (GIF.etaBin*GIF.EtaPattern)';
tempGamma = (GIF.threshold.gammaBin*GIF.threshold.GammaPattern)';
tempI = I(1,:);
tempSamplingFreq = samplingFreq;

% SIMULATE POOL
%     tempNbrRepet = nbrRepetPerPool;
%     nuGIF = nan(nbrPool,length(tempI));
%     matlabpool(nbrPool)
%     parfor p=1:nbrPool
%         nuGIF(p,:) = IFEtaMTNu(tempParam,tempEta,tempGamma,tempI,tempNbrRepet,tempSamplingFreq);
%     end
%     matlabpool close
%     nuGIF = mean(nuGIF,1);
% CLOSE POOL
nuGIF = IFEtaMTNu(tempParam,tempEta,tempGamma,tempI,nbrRepet,tempSamplingFreq);

mGIF = inprodGamma(nuGIF,nuGIF,delta);
MdStarGIF = (2*inprodGamma(nuData,nuGIF,delta)) / (mData + mGIF);
sprintf('LIF+eta+gamma => Md* = %f',MdStarGIF)

GIF.threshold.MdStarGIF = MdStarGIF; GIF.MdStarGIF = MdStarGIF;
GIF.deltaMd = delta*dt;
GIF.CPUTimeMd = toc;
GIF.nbrRepetTotal = nbrRepet;
    
fileName = sprintf('GIFTraining_%s',cName);
save(fileName,'GIF','GIFRef')

end