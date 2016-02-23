function [gradtheta Htheta] = dLdtheta(X,theta,spike,dt)
% compute likelihood of the data theta, spike

lambda0 = 1e-3;

ind = spike==1;
indNo = spike==0;

XSpike = X(:,ind);
XNoSpike = X(:,indNo);
expNoSpike = exp(theta'*XNoSpike);
clear ind indNo spike X theta

% 1. GRADIANT

gradtheta = sum(XSpike,2) - lambda0*dt*(expNoSpike*XNoSpike')';

% 2. HESSIAN

temp = nan(size(XNoSpike));
for j=1:size(XNoSpike,1)
    temp(j,:) = XNoSpike(j,:).*expNoSpike;
end
Htheta = -lambda0*dt*(temp*XNoSpike'); clear temp

clear XSpike XNoSpike expNoSpike

end