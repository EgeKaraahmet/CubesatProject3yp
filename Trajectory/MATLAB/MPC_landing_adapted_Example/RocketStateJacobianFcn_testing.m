function [A, B] = CubeSatStateJacobianFcn(x,u)
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
GM = mi*1e9;
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

%% Jacobian 
df3dh = 0.5 * rho * (Cd/m/SH) * u - 2 * grav * sin(gamma) / (Re+h); 
df3dV = -rho * (Cd/m) * V; 
df3dgamma = grav * cos(gamma); 
df3du = -0.5*rho*V^2*Cd/m; 

df4dh = -V*cos(gamma)/((Re+h)^2);
df4dV = cos(gamma)/(Re+h);
df4dgamma = -V*sin(gamma)/(Re+h);

df5dh = cos(gamma) * ((V^2*(h+Re)-2*GM)/V*(h+Re)^3);
df5dV = -cos(gamma) * (Gm + V^2*(h+Re))/(V^2*(h+Re)^2);
df5dgamma = -grav/V*sin(gamma) + V*sin(gamma)/(Re+h);



A = [   0      0  -sin(gamma)   0        0        ;
        0      0  -cos(gamma)   0        0        ;
      df3dh    0    df3dV       0    df3dgamma    ;
      df4dh    0    df4dV       0    df4dgamma    ;
      df5dh    0    df5dV       0    df5dgamma    ]; 

B = [ 0;
      0;
     df3du;
      0;
      0];