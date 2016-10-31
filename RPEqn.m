%**************************************************************************
%
% Rayleigh-Plesset Equation
%
%   Function returns the Rayleigh-Plesset equation (to be solved with an
%   ODE solver), for a shelled bubble containing an ideal gas. Includes
%   effects of surface tension and shell behavior, but not yet the buckling
%   described in the paper.
%
%   The expression is given by Eq. (3) of "A model for large amplitude 
%   oscillations of coated bubbles accounting for buckling and rupture",
%   Marmottant et al., JASA 118(6) 2005.
%
% Inputs
%   t - Time vector [s]
%   y - 
%   f_excit - Excitation frequency (sinusoidal) [Hz]
%   R0 - Equilibrium bubble radius [m]
%   pAmp - Excitation amplitude [Pa]
%
% Returns
%   dy - Vector of ODEs for solver
%     dy(1) = Rdot
%     dy(2) = Rddot
%
%**************************************************************************

function dy = RPEqn(t, y, medium, bubble, sim, excitation)

% Interpolate the excitation signal to the evaluation times
pA = interp1( excitation.tVector, excitation.signal, t );

R = y(1);
dR = y(2);

% Get medium properties
k = medium.k; % Adiabatic constant
sigma = medium.sigma; % Surface tension
rho = medium.rho; % Density
c0 = medium.c0; % Sound speed
mu = medium.mu; % Viscosity
p0 = medium.p0;

% Get bubble properties
Pv = bubble.Pvap;
R0 = bubble.R0;

% Shorthand
P_Ge = 2*sigma/R0 + p0 - Pv;

ddR = (1./R).*( ...
       (1./rho).*( ...
          P_Ge*(R/R0).^(-3*k).*(1 - 3*k*dR./c0) ...
          - p0 - 2*sigma./R - 4*mu.*dR./R - pA - Pv ...
       ) ...
       - 3/2*dR^2 ...
    );

dy = [dR; ddR];

end