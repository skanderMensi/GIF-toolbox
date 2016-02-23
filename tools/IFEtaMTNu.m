function [nu vv vThrehold] = IFEtaMTNu(param,eta,gamma,I,nbrRepet,samplingFreq)

dt = 1e3/samplingFreq;
lambda0 = 1/samplingFreq;

C = param(1); gL = param(2);
EL = param(3); Er = param(4);
tRefr = param(5)/dt;
v0 = param(6); DV = param(7);
T = length(I);
nu = zeros(T,1);

for i=1:nbrRepet
%     disp(nbrRepet+1-i)
    w = zeros(T+length(eta)+1+tRefr,1);
    VT = zeros(T+length(gamma)+1+tRefr,1);
    vv = nan(T+length(eta)+1+tRefr,1);
    vThrehold = nan(T+length(eta)+1+tRefr,1);
    
    tSpike = -T;
    v=EL;
    t=1;
    while(t<=T)
        
        dv = (-gL*(v-EL) + I(t) - w(t))/C;
        v = v + dt*dv;
        vt = v0 + VT(t);

        p = exp(-lambda0*exp((v-vt)/DV));
    
        if(rand()>p && t-tSpike>tRefr)
            nu(t) = nu(t) + 1;
            vThrehold(t) = vt;
            tSpike = t;
            
            vv(t) = v;
            vv(t+1) = 30;
            vv(t+2:t+tRefr) = Er;
            
            t = t+tRefr;
            v = Er;
            
            vv(t) = v;

            t=t+1;
            w(t:t+length(eta)-1) = w(t:t+length(eta)-1) + eta;
            VT(t:t+length(gamma)-1) = VT(t:t+length(gamma)-1) + gamma;
        else
            vv(t) = v;
            
            t=t+1;
        end
    end
end
nu = nu'/nbrRepet;

nu = nu(1:length(I));
vv = vv(1:length(I));
vThrehold = vThrehold(1:length(I));

end