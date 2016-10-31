clear all
close all
clc


% Set medium properties
medium.mode = 'static';
medium.p0 = 1.0e5;
medium.pC = 1e5;
medium.c0 = 1485;
medium.rho = 998;
medium.k = 1.33;
medium.sigma = 0.0725;  % @ 25C
medium.mu = 0.001; % @ 25C

% Set bubble properties
bubble.R0 = 2.0E-6; % Equilibrium radius
bubble.dR0 = 0; % Initial wall velocity [m/s]
bubble.Pvap = 2.33e3; %

% Set simulation properties
sim.tMin = 0;
sim.tMax = 5E-4;
sim.df = 1/sim.tMax*2;
sim.tMinMax = [sim.tMin sim.tMax];
sim.initial = [1*bubble.R0; bubble.dR0];
sim.initial0 = [.00*bubble.R0; bubble.dR0];
sim.model = 'Rayleigh-Plesset';

% Natural frequency for Rayleigh-Plesset
sim.f0_rp = 1./(2*pi*bubble.R0).*sqrt( ...
    (1./medium.rho).*(3.*medium.k*(medium.p0 ) ...
    + 2*medium.sigma./bubble.R0.*( 3*medium.k - 1) ) ...
    );
% Minneart natural frequency
sim.f0_m = 1/(2.*pi)*sqrt( ...
    3.*medium.k*medium.p0./(medium.rho.*bubble.R0.^(2)) ...
    );

sim.fo_L = 1./(2*pi*bubble.R0*sqrt(medium.rho))*sqrt( ...
    3.*medium.k*(medium.p0 + 2*medium.sigma./bubble.R0 - bubble.Pvap) ...
    - 2*medium.sigma./bubble.R0 ...
    - 4*medium.mu^2/( medium.rho*bubble.R0.^(2) ) ...
    );
%% Excitation function properties
pAmp = 50E3; % [Pa]
f0 = sim.f0_rp; % [Hz]
omega0 = 2.*pi.*f0;
tSim = linspace( sim.tMin, sim.tMax, 1000);

excitation.signal = pAmp.*( ...
    sin( omega0.*tSim ) + sin( 2.*omega0.*tSim + pi./10 ) ...
    );
excitation.tVector = tSim;

%% Main computation block

% Solve to Rayleigh-Plesset Equation
tSpan = [ sim.tMin, sim.tMax ];
initialConditions = sim.initial';

[time, solution] = ode15s(  ...
    @(t, y) RPEqn(t, y, medium, bubble, sim, excitation), ...
    tSpan, initialConditions );

%% Plot

RNorm = solution(:, 1)./bubble.R0;
t_ms = time.*1E3;

figure()
subplot( 2, 1, 1);
s = excitation.signal;
plot( excitation.tVector.*1E3, s./max(abs(s)), 'k' );
ylabel( '$s(t)$', 'FontSize', 22, 'Interpreter', 'latex' );
set( gca, 'XTickLabel', '' );

subplot( 2, 1, 2);
plot(t_ms, RNorm, 'k');

ylabel('$R/R_{0}$', 'FontSize', 22, 'Interpreter', 'latex');
xlabel('Time [ms]', 'FontSize', 18);
