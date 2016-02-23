function [spike shiftDt] = ExtractSpikeTime(voltage,VLimit,samplingFreq)
%
%   Given a voltage trace 'voltage', detect the spiketimes with the upward
%   zero crossing cirterion
%   shiftDt is the shift in timestep so that the spikes occur at the
%   maximum of the mean third drivative (in a 1 ms windows before zero crossing).
%
%   Output: spike is a vector containing 1 when a spike is present
%           shiftDt the shift in timestep to apply to spike for the third
%           derivative criterion

dt = 1e3/samplingFreq;
diffV = [diff(voltage) 0];

%VLimit = -20;                 %zero-crossing criterion
Dt = round(1/dt);           %1ms window
tRefr = floor(4.2/dt);        %refractory period of 4 ms
spike = zeros(length(voltage),1);

pattern1 = zeros(1,Dt+1);
t=1; k=0;                   %number of spikes
while(t<=length(voltage))
    if(voltage(t) >= VLimit && diffV(t) > 0)
        k = k+1;
        spike(t) = 1;
        pattern1 = pattern1 + voltage(t-Dt:t);
        t = t + tRefr - 1;
    end
    t=t+1;
end

pattern1 = pattern1/k;
pattern1 = diff(diff(diff(pattern1)));
[temp shiftDt] = max(pattern1); clear temp
shiftDt = Dt-shiftDt;

end