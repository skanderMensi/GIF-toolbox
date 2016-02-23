function [I] = noisySinWave(mu,sigma,dSigma,f,tau,T,samplingFreq)

dt = 1/samplingFreq;
Tdt = round(T*samplingFreq);
tau = tau/1e3; % correlation time constant in s

% OU NOISE %
ouNoise = nan(1,Tdt);
ouNoise(1) = 0; % initialize OU process 1
noise1 = randn(Tdt,1);

for t=1:1:Tdt-1 % forward Euler to integrate white noise and produces correlated noise
    sigmaTemp = sigma*(1 + dSigma*sin(2*pi*f*(t/samplingFreq)));
    ouNoise(t+1) = ouNoise(t) - dt*(ouNoise(t)/tau) + sigmaTemp*sqrt(2*dt/tau)*noise1(t);
end

I = mu + ouNoise;

end