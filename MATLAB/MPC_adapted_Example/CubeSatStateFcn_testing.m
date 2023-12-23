function dxdt = CubeSatStateFcn_testing(x,u)
%
% x: (1) h: altitude 
%    (2) X: distance travelled
%    (3) V: velocity 
%    (4) theta: true anomaly
%    (5) gamma: flight path angle
%
% u: fin areas
%
% The continuous-time model is valid only if the rocket above or at the
% ground (y>=10).


%constants
mi=398600.44;                  %km^3/s^2   % Earth G*M
m = 6;                         % kg
Re=6371.0088 * 10^3;           %m         % Earth mean radius 
g0=9.80665;                    %m/s^2      % Gravitational acceleration (sea level)
Cd = 2.2;

%% State 
h = x(1);
X = x(2); 
V = x(3);
theta = x(4);
gamma = x(5); 

%% Input 
% u; 

%% density 
SH = 8397.5; 
rho = 1.225 * exp(- h/(SH));    % my Python model, using SH = 8.43 and rho0=1.221 as constant
    
%A warning is printed if a negative altitude is predicted by %the simulation (due to Simulink discrete stopping criterion)
if h<0
   % fprintf('Warning: mismatched density, altitude: %3.2e km',h) 
   rho=1.225;
end

%% Gravity acceleration 
grav = mi*1e9 / (Re+h)^2;

%% d(state)/dt
doth = -V*sin(gamma);
dotX = -V*cos(gamma);
dotV = -0.5*rho*V^2*Cd*u/m + grav * sin(gamma);
dottheta = V * cos(gamma) / (Re+h); 
dotgamma = grav/V*cos(gamma) - dottheta; 


dxdt = [doth;
        dotX;
        dotV;
        dottheta;
        dotgamma;];


