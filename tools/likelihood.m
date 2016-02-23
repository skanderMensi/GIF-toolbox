function [L] = likelihood(V,Vt,spike,DeltaV,dt)
% compute likelihood of the data Vt, spike

lambda0 = 1e-3;

noSpike = 1-spike;

Lspike = (1/DeltaV)*(spike'*(V-Vt));
LnoSpike = lambda0*dt*noSpike'*exp((V-Vt)/DeltaV);

L = Lspike - LnoSpike;

end