%% 
clear 
clc
close all

%% Initialisation 
% Cubesat
M = 3   ;  % mass of the Cubesat
S = 1   ;   % Total surface area! need to consider aerobrake

% Data of the Earth
R_Earth    = 6371e3;            % Radius of the Earth (km)
GM         = 3.986004418e14;    % m^3/s^(2)
omega = [0;0;7.2921159e-5];     % Earth's rotation rate in rad/s

% Retrive from class
earth      = Earth(R_Earth, GM); 
database   = MSISE_90(); % solar activity counted data


% Speed of sound 
a = 343;            % m/s 

% Cb = Cd * S / M;  % ! This is not correct 
Cb = 200; 


% Initial states 
r0 = 120e3;         % initial radial distance (m)
v0 = 7200;          % initial velocity (m/s)



% Initial latitude and longitude
lat0 = asin(0);     % Cubesat starts from the Equator 
lon0 = 0;           % Cubesat starts from the Prime Meridian

% Convert polar coordinates to Caresian components 
[x0, y0, z0] = earth.cartesian(lat0,lon0,r0); 
[vx0,vy0,vz0] = earth.cartesian(lat0,lon0,v0);
[ax0,ay0,az0] = earth.cartesian(0,0,0); 



max_it = 1000000;  % Max iteration
h = 0.01;          % step size


% Size vector 
o = zeros(1,max_it);
x = o; y = o; z = o;
vx = o; vy = o; vz = o;
ax = o; ay = o; az = o;
t = o;



% ODE solver (Euler's method)
% Initial values
x(1) = x0; y(1) = y0; z(1) = z0;
vx(1) = vx0; vy(1) = vy0; vz(1) = vz0;

% ODE solver (Euler's method)
% Initial values
x(1) = x0; y(1) = y0; z(1) = z0;
vx(1) = vx0; vy(1) = vy0; vz(1) = vz0;

for k = 1:(length(o)-1)
    % Position vector
    p = [x(k); y(k); z(k)];
    r = norm(p);

    % Density at a given altitude
    rho = earth.density(r);

    % Relative velocity
    v = [vx(k); vy(k); vz(k)];

    % Absolute velocity
    % V = v - cross(omega, p);

    % Acceleration
    % a_aero = 0.5 * rho * Cb * norm(V) * V;
    a_aero = 0.5 * rho * Cb * norm(v) * v;
    a_g = earth.gravity(r) * (p / r);

    a = -(a_aero + a_g); 

    % Update velocity
    v = v + a * h;  % Update velocity based on acceleration
    vx(k + 1) = v(1);
    vy(k + 1) = v(2);
    vz(k + 1) = v(3);

    % Update position
    p = p + v * h;  % Update position based on velocity
    x(k + 1) = p(1);
    y(k + 1) = p(2);
    z(k + 1) = p(3);

    if earth.altitude(p) <= 0
        fprintf('Done\n');
        break;
    end
end

% Plot altitude against time
figure;
plot((0:(k-1)) * h, (sqrt(x(1:k).^2 + y(1:k).^2 + z(1:k).^2) - earth.radius) / 1000);
xlabel('Time (s)');
ylabel('Altitude (km)');
title('Altitude vs. Time');
grid on;



% Display grid lines
grid on;














