%% 
clear 
clc
close all

%% Initialisation 

% Data of the Earth
R_Earth    = 6371e3;            % Radius of the Earth (km)
GM         = 3.986004418e14;    % m^3 s^(-2)
omega = 7.2921159e-5;           % Earth's rotation rate in rad/s

% Retrive from class
earth      = Earth(R_Earth, GM); 


% Initial states 
r0 = 120e3;         % initial radial distance (m)
v0 = 7.2;          % initial velocity (m/s)



% Initial latitude and longitude
lat0 = asin(0);     % Cubesat starts from the Equator 
lon0 = 0;           % Cubesat starts from the Prime Meridian
% Convert polar coordinates to Caresian components 
[x0, y0, z0] = earth.cartesian(lat0,lon0,r0); 
[vx0,vy0,vz0] = earth.cartesian(lat0,lon0,v0);
[ax0,ay0,az0] = earth.cartesian(0,0,0); 


% Time span
tspan = [0, 1000];  % Time span in seconds

% Initial state vector
initial_state = [x0, y0, z0, vx0, vy0, vz0];

%% Solver 
% Solve the system of equations
[t, state] = ode45(@motion_eqns, tspan, initial_state);


%% Results 
% Extract the state variables
x = state(:, 1);
y = state(:, 2);
z = state(:, 3);
vx = state(:, 4);
vy = state(:, 5);
vz = state(:, 6);

%%

% Initialize arrays to store altitude and time
altitude = zeros(size(state, 1), 1); % Initialize with zeros
time = t;

% Loop through each time step to convert and calculate altitude
for i = 1:size(state, 1)
    % Extract position and velocity components
    x = state(i, 1);
    y = state(i, 2);
    z = state(i, 3);
    vx = state(i, 4);
    vy = state(i, 5);
    vz = state(i, 6);

    % Convert position and velocity to polar coordinates
    [lat, lon, alt] = earth.polar(x, y, z);
    [vlat, vlon, valt] = earth.polar(vx, vy, vz);

    % Calculate altitude from spherical coordinates
    altitude(i) = alt;

end

% Plot altitude against time
figure;
plot(time, altitude);
xlabel('Time (s)');
ylabel('Altitude (m)');
title('Altitude vs. Time');
grid on;





    