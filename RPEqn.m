%**************************************************************************
%
% Rayleigh-Plesset Equation (Marmottant Model)
%
%   Function returns the state vector (to be passed to an ODE solver), for
%   a shelled bubble containing an ideal gas. The surface tension is may be
%   a function of radius, for an arbitrary input excitation pulse.
%
%   Based on the modified Rayleigh-Plesset equation given in Ref. [1].
%
% Inputs
%   t - Time vector [s]
%   y - State vector
%   medium. Struct describing the (homogeneous) medium:
%     p0    - Ambient pressure [Pa]
%     c0    - Sound speed [m/s]
%     rho   - Density [kg/m^3]
%     k     - Polytropic exponent
%     sigma - Surface Tension [N/m]
%     mu    - Dynamic viscosity [Pa s]
%   bubble. - Struct describing bubble parameters
%     R0       - Equilibrium radius [m]
%     Rbuckle  - Buckling radius [m]
%     Rbreak   - Breaking radius [m]
%     Pvap     - Vapor pressure of contained gas [Pa]
%     hasShell - 0 if bubble is free, 1 if it has a shell
%     shell.   - If bubble.hasShell ~= 0:
%       thickness     - Shell thickness [m]
%       bulkViscosity - Bulk viscocity of shell material [Pa s]
%   excitation.  Struct containing excitation vectors
%     tVector - Time vector [s]
%     signal  - Pressure incident on the bubble [Pa]
%
% Returns
%   dy - Vector of ODEs for solver
%     dy(1) = Rdot
%     dy(2) = Rddot
%
%**************************************************************************

function dy = RPEqn(t, y, medium, bubble, excitation)

% Interpolate the excitation signal to the evaluation times t
pA = interp1( excitation.tVector, excitation.signal, t );

% Get variables from state vector
R  = y(1);
dR = y(2);

% Get medium properties
k = medium.k;         % Adiabatic constant
sigma = medium.sigma; % Surface tension
rho = medium.rho;     % Density
c0 = medium.c0;       % Sound speed
mu = medium.mu;       % Viscosity
p0 = medium.p0;

% Get bubble properties
Pv = bubble.Pvap;
R0 = bubble.R0;

% Set shell properties if they exist
if bubble.hasShell
    epsilon = bubble.shell.thickness;
    muShell = bubble.shell.bulkViscosity;
    kappa_s = 3.*epsilon.*muShell;
else
    kappa_s = 0;
end

% Compute expression for second derivative of the bubble radius [Eq. (3)
% in Ref. 1].
ddR = (1./R)*( ...
       (1./rho)*( ...
          (p0 - Pv + 2*sigma(R0)/R )*(R/R0)^(-3*k)*(1 - 3*k*dR/c0) ...
          - p0 - 2*sigma(R)/R - 4*mu*dR/R - 4*kappa_s*dR/R^(2) - pA - Pv ...
       ) ...
       - 3/2*dR^(2) ...
    );

% Return the state vector
dy = [dR; ddR];

end

% Citations
% 1. Marmottant et. al., "A model for large amplitude oscillations of 
%      coated bubbles accounting for buckling and rupture", J. Acous. Soc.
%      Am. 118(6) pp. 3499--3505 (2005).