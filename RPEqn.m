function dy = RPEqn(t, y, f_exc,R0,pAmp)

global medium simul bubble

R = y(1);
dR = y(2);

switch medium.mode
    case 'heaviside'
        t0 = 500e-6;
        if t<t0
            p0 = medium.p0;
        elseif t>=t0
            
            p0 = medium.p0 + medium.pC;
        end
        
    case 'static'
        p0 = medium.p0;
        
end

k = medium.k;
sigma = medium.sigma;
rho = medium.rho;
c0 = medium.c0;
mu_ = medium.mu;
Pv = bubble.Pvap;

P_Ge = 2*sigma/R0 + p0 - Pv;
pA_ = simul.excite.polarity*pAmp*sin(2*pi*f_exc*t);
if strcmp(simul.excite.type,'Pulsed')
    if t<= 1/f_exc*simul.excite.numcycles
        pA_ = simul.excite.polarity*pAmp*sin(2*pi*f_exc*t);
    else
        pA_ = 0;
    end
end

ddR = ((P_Ge*(R0/R)^(3*k)*(1) - 2*sigma/R - 4*mu_*dR/R - p0 - pA_ + Pv)/rho - 3/2*dR^2)/R;
% ddR = ((P_Ge*(R0/R)^(3*k)*(1-3*k*dR/c0) - 2*sigma/R - 4*mu_*dR/R - p0 - pA_ - Pv)/rho - 3/2*dR^2)/R;

dy = [dR; ddR];