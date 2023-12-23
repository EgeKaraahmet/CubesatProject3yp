function [A, B] = CubeSatStateJacobianFcn_23122023(x,u)
% x: (1) X: x-component of position vector p 
%    (2) Y: y-component of position vector p 
%    (3) Vx: x-component of velocity vector V
%    (4) Vy: y-component of velocity vector V
%
% u: fin areas
%
% The continuous-time model is valid only if the rocket above or at the
% ground (y>=10).


%constants
m = 6;                         % kg
GM = 398600.44*1e9;            %m^3/s^2   % Earth G*M
Re=6371.0088 * 10^3;           %m         % Earth mean radius 
g0=9.80665;                    %m/s^2      % Gravitational acceleration (sea level)
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
df3dX = (C1*Cd*V1*v_mag*X*rho_exp)/(u*C2*m*r) + 3*GM*X^2/r^5 - GM/r^3 ; 
df3dY = (C1*Cd*V1*v_mag*Y*rho_exp)/(u*C2*m*r) + 3*GM*X*Y/r^5 ; 
df3dV1 = -C1*Cd*V1.^2*rho_exp/(u*m*v_mag) - C1*Cd*v_mag*rho_exp/(u*m);
df3dV2 = -C1*Cd*V1*V2*rho_exp/(u*m*v_mag); 
df3du = C1*Cd*V1*v_mag*rho_exp/(u.^2*m); 


df4dX = (C1*Cd*V2*v_mag*X*rho_exp)/(u*C2*m*r) + 3*GM*X*Y/r^5 ; 
df4dY = (C1*Cd*V2*v_mag*Y*rho_exp)/(u*C2*m*r) + 3*GM*Y^2/r^5 - GM/r^3 ; 
df4dV1 = -C1*Cd*V1*V2*rho_exp/(u*m*v_mag); 
df4dV2 = -C1*Cd*V2.^2*rho_exp/(u*m*v_mag) - C1*Cd*v_mag*rho_exp/(u*m);
df4du = C1*Cd*V2*v_mag*rho_exp/(u.^2*m); 


A = [0       0     1      0   ;
     0       0     0      1   ;
     df3dX df3dY df3dV1 df3dV2;
     df4dX df4dY df4dV1 df4dV2];

B = [0;
     0;
     df3du;
     df4du]; 




