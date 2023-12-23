function dxdt = CubeSatStateFcn_23122023(x,u)
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
m = 6;                         % kg
GM = 398600.44*1e9;            %m^3/s^2   % Earth G*M
Re = 6371.0088 * 10^3;           %m         % Earth mean radius 
g0 = 9.80665;                    %m/s^2      % Gravitational acceleration (sea level)
Cd = 2.2;

%% State 
X = x(1);
Y = x(2);
V1 = x(3); 
V2 = x(4); 

% For simplification purposes 
r = sqrt(X.^2 + Y.^2);
h = r - Re; 
v_mag = sqrt(V1.^2 + V2.^2); 



%% Input 
% u; 
beta = m / (Cd * u); 

%% density 
C1 = 0.5 * 1.225; 
C2 = 8397.5; 
% rho = 1.225 * exp(- h/(C2));    % my Python model, using SH = 8.43 and rho0=1.221 as constant    
rho_exp = exp((-h) / C2);  



%% d(state)/dt
dotX = V1;
dotY = V2;
dotV1 = -GM*X/(r^3) - C1 * rho_exp * v_mag * V1 / beta; 
dotV2 = -GM*Y/(r^3) - C1 * rho_exp * v_mag * V2 / beta; 

dxdt = [dotX;
        dotY;
        dotV1;
        dotV2];





