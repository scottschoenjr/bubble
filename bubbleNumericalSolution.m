global medium simul bubble

%% Medium properties
medium.mode = 'static';
medium.p0 = 1.0e5;
medium.pC = 1e5;
medium.c0 = 1485;
medium.rho = 998;
medium.k = 1.33;
medium.sigma = 0.0725;  % @ 25C
medium.mu = 0.001; % @ 25C

%% Medium and bubble properties
bubble.R0 = 2.0e-6;
bubble.dR0 = 0;
bubble.Pvap = 2.33e3;

%% Simulation properties
simul.tMin = 0;
simul.tMax = 1e-3;
simul.df = 1/simul.tMax*2;
simul.tMinMax = [simul.tMin simul.tMax];
simul.initial = [1*bubble.R0; bubble.dR0];
simul.initial0 = [.00*bubble.R0; bubble.dR0];
simul.model = 'Rayleigh-Plesset';

%% Natural frequency for Rayleigh-Plesset
simul.fo_rp = 1./(2*pi*bubble.R0).*sqrt(1/medium.rho*(3*medium.k*(medium.p0 )...
    + 2*medium.sigma./bubble.R0*(3*medium.k-1)));
simul.fo_l = 1/2/pi*sqrt(3*medium.k*medium.p0/medium.rho/bubble.R0^2);

simul.fo_L = 1./(2*pi*bubble.R0*sqrt(medium.rho))*sqrt(3*medium.k*(medium.p0+2*medium.sigma./bubble.R0-bubble.Pvap)...
       - 2*medium.sigma./bubble.R0-4*medium.mu^2/(medium.rho*bubble.R0^2));
%% Excitation function properties
% simul.excite.f = 20.0*1e3; % (Hz)
simul.excite.f =0.15*simul.fo_L:simul.fo_L/20:2.5*simul.fo_L; %(Hz)
simul.excite.pA = 50*1e3;  % (Pa)
excitation_type_avail = {'CW', 'Pulsed'};
simul.excite.type = excitation_type_avail{1};
simul.excite.polarity = 1;

% number of cycles
if strcmp(simul.excite.type,'Pulsed')
    simul.excite.numcyc = 4;
end

%% Main computation block
s1 = numel(simul.excite.f);
s2 = numel(simul.excite.pA);

soln = cell(s1,s2,1);

result.R_R0 = zeros(s1,s2);
result.resp = zeros(s1,s2);

f1 = figure;	ph1 = plot(simul.tMin,bubble.R0);
ah1 = gca;  set(ah1,'XLim',[simul.tMin simul.tMax]);

for i1=1:s1
    for i2=1:s2
        % Rayleigh-Plesset Equation
        a = ode15s(@(t,y)RPEqn(t, y, simul.excite.f(i1), bubble.R0,simul.excite.pA(i2)), ...
            [simul.tMin simul.tMax], simul.initial');
        result.resp(i1,i2)=(max(a.y(1,end-1000:end))-bubble.R0)/bubble.R0;
        % plot data
        figure(1)
        plot( simul.excite.f(1:i1)/simul.fo_L, result.resp(1:i1),'*')
        title(['P_A = ',num2str(simul.excite.pA(i2)/1e3),' kPa  &  F_exc = ',num2str(simul.excite.f(i1)/1e6),' MHz']),drawnow
    end
end