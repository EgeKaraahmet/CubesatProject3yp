function xkplus1 = CubeSatStateFcn_discrete_testing(x,u)
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
hk = x(1);
Xk = x(2); 
Vk = x(3);
thetak = x(4);
gammak = x(5); 

%% Input 
% u; 
uk = u; 

%% density 
SH = 8397.5; 
rho = 1.225 * exp(- hk/(SH));    % my Python model, using SH = 8.43 and rho0=1.221 as constant
    
%A warning is printed if a negative altitude is predicted by %the simulation (due to Simulink discrete stopping criterion)
if hk<0
   % fprintf('Warning: mismatched density, altitude: %3.2e km',h) 
   rho=1.225;
end

%% Gravity acceleration 
grav = mi*1e9 / (Re+hk)^2;


%% Discretisation (Euler's forward method)
step = 0.01; 

hkplus1 = hk + step*(-Vk*sin(gammak));
Xkplus1 = Xk + step*(-Vk*cos(gammak));
Vkplus1 = Vk + step*(-0.5*rho*Vk^2*Cd*uk/m + grav * sin(gammak));
thetakplus1 = thetak + step*(Vk * cos(gammak) / (Re+hk)); 
gammakplus1 = gammak + step*(grav/Vk*cos(gammak) - thetakplus1); 


xkplus1 = [hkplus1;
        Xkplus1;
        Vkplus1;
        thetakplus1;
        gammakplus1;];


