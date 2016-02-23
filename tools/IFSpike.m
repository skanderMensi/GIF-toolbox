function [V] = IFSpike(param,I,spike,samplingFreq)

dt = 1e3/samplingFreq;
C = param(1); gL = param(2);
EL = param(3); Er = param(4);
tRefr = param(5)/dt;
T = length(I);
V = nan(T,1);

spike = [spike+1;inf];
v=EL;
ind = 1;
t=1;
while(t<=T)
    dv = (-gL*(v-EL) + I(t))/C;
    v = v + dt*dv;
    
    V(t) = v;
    
    if(t==spike(ind))
        V(t) = v;
        V(t+tRefr) = Er;
        t = t+tRefr;
        v = Er;
        ind = ind + 1;
    end
    t=t+1;
end
V = V(1:length(I));
end