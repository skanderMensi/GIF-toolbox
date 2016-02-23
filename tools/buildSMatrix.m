function [S SPattern tBin] = buildSMatrix(nbrBin,T,spike,tRefr,kernelType,samplingFreq)

dt = 1e3/samplingFreq;

if(strcmp(kernelType,'eta'))
    bin0 = 2/dt;
elseif(strcmp(kernelType,'gamma'))
%     bin0 = 10/dt;
    bin0 = 2/dt;
elseif(strcmp(kernelType,'gammaFS'))
    bin0 = 5/dt;
end

binSize = ceil(diff(logspace(log10(1), log10(T/dt), nbrBin+1)))+bin0;
binSize(end) = binSize(end) - abs((T/dt)-sum(binSize));

SPattern = zeros(nbrBin,round(T/dt));
tBin = nan(nbrBin,1);
for i=1:nbrBin
    if(i==1)
        SPattern(i,1:binSize(1)) = 1;
        tBin(i) = round(binSize(i)/2)*dt;
    else
        SPattern(i,sum(binSize(1:i-1))+1:sum(binSize(1:i))) = 1;
        tBin(i) = (sum(binSize(1:i-1))+1+round(binSize(i)/2))*dt;
    end
end

S = zeros(nbrBin,length(spike)+round(T/dt)+tRefr/dt);
spiketime = (find(spike==1));
spiketime = spiketime + round((tRefr)/dt);
for i=1:length(spiketime)
    S(:,spiketime(i)+1:spiketime(i)+round(T/dt)) = ...
        S(:,spiketime(i)+1:spiketime(i)+round(T/dt)) + SPattern;
end

S = S(:,1:length(spike));

end