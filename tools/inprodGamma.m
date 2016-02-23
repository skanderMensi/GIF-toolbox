function [output] = inprodGamma(x,y,delta)
%
%   Compute inner product Gamma with replacement according to
%   <X,Y> = int_{0}^{T} int_{-inf}^{inf} int_{-inf}^{inf} ...
%   K_{Delta}(s,s')X(t-s)Y(t-s')dsds'dt
%   with K_{Delta}(s,s') = the coincidence detector windows (i.e a rect)
%
%   Here we have:
%   K_{Delta}(s,s') = dirac(s)h_{r}(s';Delta)
%   h_{r}(s;Delta) = Heaviside(s-Delta)Heaviside(Delta-s)
%   
%   X and Y are N-dimensionnal vectors

filt = ones((2*delta)+1,1);

x1 = fftfilt(filt,[x zeros(1,length(filt))]);
x1 = x1(delta+1:end); x1 = x1(1:length(x));
output = sum(x1.*y);

end