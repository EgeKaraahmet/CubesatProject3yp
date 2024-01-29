function dxdt = CubeSatStateFcn_25122023(x,u)
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
mi=398600.44;                  % m^3/s^2   % Earth G*M
m = 6;                         % kg
Re=6371.0088 * 10^3;           % m         % Earth mean radius 
g0=9.80665;                    % m/s^2      % Gravitational acceleration (sea level)
Cd = 2.2;
%space environment      
date=[2023,01,01,12,00,00];         % intial time of simulation [y,m,d,h,m,s]
jdate=juliandate(date);
F107_avg=90.85;         %SFU        % F10.7 average of 3x27 days before the date under consideration
F107_day=86.0;          %SFU        % F10.7 average of day before the date under consideration


Kp=1;                                 % Kp three-hourly planetary geomagneticindex

%% State 
h = x(1);         % m 
X = x(2);         % m 
V = x(3);         % m/s
theta = x(4);
gamma = x(5); 

%% Input 
% u; 

%% density 
% 
SH = 8397.5; 
rho = 1.225 * exp(- h/(SH));    % my Python model, using SH = 8.43 and rho0=1.221 as constant

%A warning is printed if a negative altitude is predicted by %the simulation (due to Simulink discrete stopping criterion)
if h<0
   % fprintf('Warning: mismatched density, altitude: %3.2e km',h) 
   rho=1.225;
end

% 29012024 updates
% [~,~,rho,~,~] = ATM1976(h);



%% Gravity acceleration 
grav = mi*1e9 / (Re+h)^2;

%% Ballistic coefficient
BC  = m / Cd / u; 

%% d(state)/dt
doth = -V*sin(gamma);
dotX = -V*cos(gamma);
dotV = -0.5*rho*V^2/BC + grav * sin(gamma);
dottheta = V * cos(gamma) / (Re+h); 
dotgamma = grav/V*cos(gamma) - dottheta; 


dxdt = [doth;
        dotX;
        dotV;
        dottheta;
        dotgamma];


