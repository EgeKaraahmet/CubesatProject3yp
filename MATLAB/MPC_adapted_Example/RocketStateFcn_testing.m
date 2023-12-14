function dxdt = RocketStateFcn(x,u)
% In a 2D environment with standard XY axis, the robot is a circular disc
% (20 meters in diamater).  Two thrusts are to the left and right of the
% center.  Tilting (theta) is defined as positive to left and negative to
% the right (0 means robot is vertical).
%
% x: (1) h: altitude 
%    (2) X: distance travelled
%    (3) V: velocity 
%    (4) Vx: horizontal velocity 
%    (5) Vy: verti velocity 
%    (6) gamma: flight path angle
%    (7) theta: true anomaly 
%
% u: fin areas
%
% The continuous-time model is valid only if the rocket above or at the
% ground (y>=10).


m = 6;     % mass of CubeSat (kg)   
g = 9.81;  % gravity (m/s^2)

% inertia for a flat disk
I = 0.5*m*L1^2;
% get force and torgue
Tfwd   = u(2) + u(1);
Ttwist = u(2) - u(1);
Fx = -sin(x(3)) * Tfwd;
Fy =  cos(x(3)) * Tfwd;
Mz =  Ttwist * L2  ;
% treat as "falling" mass
dxdt = [x(4); x(5); x(6); Fx/m; Fy/m - g; Mz/I];


