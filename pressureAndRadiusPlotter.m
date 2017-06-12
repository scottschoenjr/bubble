clear all
close all
clc

% Set bubble properties
bubble.R0 = 2.8E-6;   % Equilibrium radius [m]
bubble.dR0 = 0;       % Initial wall velocity [m/s]
bubble.Pvap = 2.33E3; % Vapor Pressure
bubble.hasShell = 0;
bubble.shell.thickness = 10E-9; % [m] Estimation!!!
bubble.shell.bulkViscosity = 1; % Bulk viscosity of shell [Pa]
bubble.Rbuckle = 0.99.*bubble.R0; % Radius at which shell buckles [m]
bubble.Rbreak = 1.01*bubble.R0; % Radius at which shell separates [m]

% Set medium properties (from Tu et. al.)
medium.p0 = 1.013E5;
medium.c0 = 1500;
medium.rho = 100;
medium.k = 1.07; 
medium.mu = 0.002;

% Define sigma as a function of R
medium.sigma = @(R) ...
    (R <= bubble.Rbuckle ).*...
      ( 0 ) + ...
    ( (R > bubble.Rbuckle) && (R <= bubble.Rbreak) ).* ...
      ( 0.55.*( (R./bubble.Rbuckle).^(2) - 1 ) ) + ...
    (R > bubble.Rbreak ).* ...
      ( 0.0725 ); 
  
% Set to constant if the bubble doesn't have a shell
if ~( bubble.hasShell )
   medium.sigma = @(R) 0.0725;    
end

% Set simulation properties
sim.tMin = 0;
sim.tMax = 2E-4;
sim.df = 1/sim.tMax*2;
sim.tMinMax = [sim.tMin sim.tMax];
sim.initial = [1*bubble.R0; bubble.dR0];

% Minneart natural frequency
sim.f0_m = 1/(2.*pi)*sqrt( ...
    3.*medium.k*medium.p0./(medium.rho.*bubble.R0.^(2)) ...
    );

%% Excitation function properties
pAmp = 1E5; % [Pa]
f0 = 0.9.*sim.f0_m; % [Hz]
omega0 = 2.*pi.*f0;
tSim = linspace( sim.tMin, sim.tMax, 5E4);
t0 = 2E-6; % Make sure not too large that ODE solver misses excitation
t1 = t0 + 6./f0;
BW = 0.2;

excitation.tVector = tSim;
excitation.signal = pAmp.*excitationPulse( tSim, f0, BW, t0, 0 );

useImpulse = 1;
if useImpulse
    s = 0.*tSim;
    startIndex = find( tSim > t0, 1 );
    endIndex = startIndex + 4;
    s( startIndex : endIndex ) = pAmp;
    excitation.signal = s;
end

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

% FFT of bubble wall velocity motion
Rtilde = fft( RNorm );
RdotTilde = fft( Rdot );
% Sampling interval is nonuniform -- Is avg best?
dt = t(2) - t(1); 
Fs = 1./dt;
fVector = linspace(0, Fs, length(t));

%% Plot

figure()
set( gcf, 'Position', [50, 50, 1000, 800] );

% Excitation
signalPlot = subplot( 2, 1, 1);
s = excitation.signal;
plot( t_us, s./max(abs(s)), 'b' );

title( sprintf( '($p_{0}$ = %3.1f kPa)', pAmp./1E3 ), 'Interpreter', 'latex'  );
xlim(1E6.*[sim.tMin, sim.tMax]);

ylabel( '$P_{\rm ac}(t)$', 'FontSize', 22, 'Interpreter', 'latex' );
set( gca, 'XTickLabel', '' );

ylim([-1.05, 1.05]);

% Radius
radiusPlot = subplot( 2, 1, 2);
plot(t_us, RNorm, 'k');

xlim(1E6.*[sim.tMin, sim.tMax]);
% ylim([0, 2]);

ylabel('$R/R_{0}$', 'FontSize', 28, 'Interpreter', 'latex');
xlabel('Time [$\mu$s]', 'FontSize', 24);

% Link axes
linkaxes( [signalPlot, radiusPlot], 'x' );
zoom xon;

% xlim([0.0, 8])

% ----- Axes Limits
% set( signalPlot, 'xlim', [8, 14] );
% set( radiusPlot, 'xlim', [8, 14], 'ylim', [0, 3.5] );

% Citations
% 1. Juan Tu et. al., "Microbubble Sizing and Shell Characterization Using
%      Flow Cytometry", IEEE Trans. Ultrason. Ferroelectr. Freq. Control 
%      58(5) pp. 955--963 (2011)
