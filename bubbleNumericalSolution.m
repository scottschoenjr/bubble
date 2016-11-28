clear all
close all
clc

% Set medium properties
medium.p0 = 1.0e5;
medium.c0 = 1485;
medium.rho = 998;
medium.k = 1.33;
medium.sigma = 0.0725;  
medium.mu = 0.001;

% Set bubble properties
bubble.R0 = 1.0E-6;   % Equilibrium radius [m]
bubble.dR0 = 0;       % Initial wall velocity [m/s]
bubble.Pvap = 2.33E3; % Vaporization Pressure
bubble.hasShell = 0;
bubble.shell.thickness = 1E-8;
bubble.shell.bulkViscosity = 1E-3; % Bulk viscosity of liquid [Pa]

% Set simulation properties
sim.tMin = 0;
sim.tMax = 2E-4;
sim.df = 1/sim.tMax*2;
sim.tMinMax = [sim.tMin sim.tMax];
sim.initial = [1*bubble.R0; bubble.dR0];

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
pAmp = 10.*101E3; % [Pa]
f0 = sim.f0_rp; % [Hz]
omega0 = 2.*pi.*f0;
tSim = linspace( sim.tMin, sim.tMax, 5E4);
t0 = 1E-5; % Make sure not too large that ODE solver misses excitation
t1 = t0 + 6./f0;

excitation.signal = pAmp.*excitationPulse( tSim, f0, 0.2, t0, 0 );
excitation.tVector = tSim;

%% Main computation block

% Solve to Rayleigh-Plesset Equation
tSpan = [ sim.tMin, sim.tMax ];
initialConditions = sim.initial';

[time, solution] = ode15s(  ...
    @(t, y) RPEqn(t, y, medium, bubble, excitation), ...
    tSpan, initialConditions );

% Get solution vectors of interets
R = solution(:, 1);    % [m]
Rdot = solution(:, 2); % [m/s]

% R and Rdot values are given at solution time points in time vector. These
% points are not uniformly spaced, so interpolate the solutions so that
% the R and Rdot values are known at the times in the excitation time vector
R = interp1( time, R, excitation.tVector );
Rdot = interp1( time, Rdot, excitation.tVector );

% Get plotting vectors
RNorm = R./bubble.R0;   % [normalized]
t = excitation.tVector; % [s]
t_us = 1E6.*t;          % [us]


% FFT of bubble motion
Rtilde = fft( RNorm );
% Sampling interval is nonuniform -- Is avg best?
dt = t(2) - t(1); 
Fs = 1./dt;
fVector = linspace(0, Fs, length(t));

%% Plot

figure()
set( gcf, 'Position', [50, 50, 1000, 800] );

% Excitation
signalPlot = subplot( 3, 1, 1);
s = excitation.signal;
plot( t_us, s./max(abs(s)), 'k' );

xlim(1E6.*[sim.tMin, sim.tMax]);

ylabel( '$s(t)$', 'FontSize', 22, 'Interpreter', 'latex' );
set( gca, 'XTickLabel', '' );

% Radius
radiusPlot = subplot( 3, 1, 2);
plot(t_us, RNorm, 'k');

xlim(1E6.*[sim.tMin, sim.tMax]);
% ylim([0, 2]);

ylabel('$R/R_{0}$', 'FontSize', 22, 'Interpreter', 'latex');
xlabel('Time [$\mu$s]', 'FontSize', 18);

% Link axes
linkaxes( [signalPlot, radiusPlot], 'x' );
zoom xon;

% ----- Axes Limits
% set( signalPlot, 'xlim', [8, 14] );
% set( radiusPlot, 'xlim', [8, 14], 'ylim', [0, 3.5] );

% Spectrum
spectrumPlot = subplot(3, 1, 3);

fNorm = fVector(2:end)./sim.f0_rp;
RtildeNorm = abs(Rtilde(2:end))./max(abs(Rtilde(2:end)));

plot( fNorm, RtildeNorm, 'k' );

xlabel( '$f/f_{0}$', 'FontSize', 22, 'Interpreter', 'latex' );
ylabel( '$\tilde{R}$', 'FontSize', 22, 'Interpreter', 'latex' );

xlim( [0, 10] );
