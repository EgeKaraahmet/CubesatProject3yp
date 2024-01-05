% ARES - Academic Re-Entry Simulator - Sept 2019
% This script evaluates the decay of an Earth orbiting spacecraft exposed
% to atmospheric drag.
% Atmospheric density is calculated through:
% -Jacchia J71 model [100 to 2500 km of altitude]
% -Exponential atmosphere [0 to 100 km of altitude]
clear
close all
clc

%% Initial conditions 
% Input 
A=125 * 10^(-4);          %m^2        % spacecraft cross - sectional area

% constant 
m = 6;                    %kg         % spacecraft mass  m = 3 for ISS
Cd=2.2;                             % drag coefficient
BC = m / Cd / A; 



%constants
mi=398600.44;           %km^3/s^2   % Earth G*M
Re=6371.0088;           %km         % Earth mean radius
g0=9.80665;             %m/s^2      % Gravitational acceleration (sea level)
%spacecraft parameters
E=0;                                % aerodynamic efficiency
CL=Cd*E;                            % lift coefficient
rcurv=0.1;              %m          % TPS nose curvature radius

%space environment      
date=[2023,01,01,12,00,00];         % intial time of simulation [y,m,d,h,m,s]
jdate=juliandate(date);
F107_avg=90.85;         %SFU        % F10.7 average of 3x27 days before the date under consideration
F107_day=86.0;          %SFU        % F10.7 average of day before the date under consideration


Kp=1;                                 % Kp three-hourly planetary geomagneticindex

alt_0 = 200;              %km         % initial altitude    %% 400 km
r_a=alt_0+Re;             %km         % radius at apoapsis
r_p=alt_0+Re;             %km         % radius at periapsis
a=(r_a+r_p)/2;            %km         % semimajor axis
e=(r_a-r_p)/(r_a+r_p);                % eccentricity


%re-entry path initial values (at t=0 s)
x0=0;                             % m          % travelled distance
gamma0=0;                         % rad        % flight path angle
theta0=0;                         % rad        % true anomaly 
r0=a*(1-e^2)/(1+e*cos(theta0));   % km         % position vector length
h0=r0-Re;                         % km         % height 
V0=sqrt(mi*(2/r0-1/a));           % km/s       % orbital speed (ellipse)
V0 = 7.788; 

%de-orbiting retrograde burn
dV=0;                             % km/s       % impulsive delta V obtained
V0=V0-dV;                         % km/s       % effective inital speed 

%integration method
solv_kep='ode4';                  % solver method for orbital phase. ode4 is runge kutta
step_kep='30';  %s                % keplerian intergartor fixed step size   % 30 
stop_h=150;     %km               % thres alt for switch from Kep. to Atm. phase 

solv_atm = 'ode4';                % solver method for re-entry phase. ode4 is runge kutta
step_atm='0.1'; %s                % atmospheric integrator fixed step size     % 0.1

%conversions to I.S.
Re=Re*1000;                     %m
V0=V0*1000;                     %m/s
h0=h0*1000;                     %m
stop_h=stop_h*1000;             %m

%% Simulation 
%a 3 DOF simulator has been implemented in Simulink model SatSim_ARES
%Kepleran phase of the simulation, slow variations in states due to low
%density, thus low drag perturbation.
%Integration step: 30 s
%First run stops when altitude reachesG stop_h threshold

kepOut = sim('reference_signal_generator','Solver',solv_kep,'FixedStep',step_kep);

%   updating initial conditions, 2nd run using ouput from 1st

%space environment
ndays1=kepOut.time(end)/60/60/24;  %days since intial time of simulation
jdate=jdate+ndays1;                %new julian date at beginning of 2nd run
F107_avg=90.85;         %SFU: updating
F107_day=86.0;          %SFU: updating
Kp=1;

%re-entry path initial values (at t=0 s)
x0=kepOut.x(end);            %m        %travelled distance
gamma0=kepOut.gamma(end);    %rad      %flight path angle
theta0=kepOut.theta(end);    %rad      %true anomaly
h0=kepOut.h(end);            %m        %height
V0=kepOut.V(end);            %m/s      %orbital speed

%Re-entry simulation, atmospheric phase, more accuracy is needed.
%Smaller integration step required (e.g. 0.1 s)
%Simulation stops when altitude reaches 0 km (default)
%%%%%%%%%%%%%%%%%%
stop_h=0;   %m          %final desired altitude
%%%%%%%%%%%%%%%%%%

atmOut = sim('reference_signal_generator','Solver',solv_atm,'FixedStep',step_atm,'StartTime','kepOut.time(end)');

%% Results analysis and plotting
%   union of 1st and 2nd run results

h_reference_signal=[kepOut.h;atmOut.h];
x_reference_signal=[kepOut.x;atmOut.x];
time_reference_signal=[kepOut.time;atmOut.time];
V_reference_signal=[kepOut.V;atmOut.V];
Vx_reference_signal=[kepOut.Vx;atmOut.Vx];
Vz_reference_signal=[kepOut.Vz;atmOut.Vz];
gamma_reference_signal=[kepOut.gamma;atmOut.gamma];
theta_reference_signal=[kepOut.theta;atmOut.theta];
rho_reference_signal=[kepOut.rho;atmOut.rho];
dotV_reference_signal=[kepOut.dotV;atmOut.dotV];

%   data filtering
%section to be revised
%excludes unsensed data due to discrete stopping criteria
%only (end) value is typically wrong (i.e. h(end)<0)


%%
n_sa=length(time_reference_signal);
k=0;
for i=n_sa:-1:1
    if h_reference_signal(i)<0 %|| h(i)>1.05*(r_a*1000-Re)
        k=k+1;
        h_reference_signal(i)=[];
        x_reference_signal(i)=[];
        time_reference_signal(i)=[];
        V_reference_signal(i)=[];
        Vz_reference_signal(i)=[];
        Vx_reference_signal(i)=[];
        gamma_reference_signal(i)=[];
        theta_reference_signal(i)=[];
        rho_reference_signal(i)=[];
        dotV_reference_signal(i)=[];
    end 
end
%%
%   data analysis
%time to deorbit
% nyears=time(end)/3.154e+7;
% fprintf('Total time to de-orbit: %3.2f years\n',nyears)
fprintf('Fin area: %3.2f cm^2\n',A*10^4)

ndays=time_reference_signal(end)/3.154e+7*365;
fprintf('Total time to de-orbit: %3.2f days\n',ndays)

%structural loading
gload=abs(dotV_reference_signal/g0);             %strutural loading [g]
gload_max=max(gload);
h_gload_max=h_reference_signal(gload==gload_max)/1000;
fprintf('Max load factor of %.1f g experienced at an altitude of %.1fkm\n',gload_max,h_gload_max)

%aerothermal load
heat1=1.83e-4*V_reference_signal.^3.*sqrt(rho_reference_signal./rcurv);
heat1_max=max(heat1);
h_heat1_max=h_reference_signal(heat1==heat1_max)/1000;
fprintf('Max heat flux of %3.2e W/m^2 experienced at an altitude of %.1fkm\n',heat1_max,h_heat1_max)
%conversions
h_reference_signal=h_reference_signal/1000;     %[m] to [km]
x_reference_signal=x_reference_signal/1000;     %[m] to [km]
time_reference_signal=time_reference_signal/60; %[s] to [min]
V_reference_signal=V_reference_signal/1000;     %[m/s] to [km/s]
Vz_reference_signal=Vz_reference_signal/1000;   %[m/s] to [km/s]
Vx_reference_signal=Vx_reference_signal/1000;   %[m/s] to [km/s]
gamma_reference_signal=rad2deg(gamma_reference_signal);  %[rad] to [deg]
%theta=rad2deg(theta);
%plotting
figure
plot(V_reference_signal,h_reference_signal)
xlabel('velocity [km/s]')
ylabel('altitude [km]')
grid on
grid minor

figure
plot(h_reference_signal,heat1)
xlabel('altitude [km]')
ylabel('heat flux [W/m^2]')
grid on
grid minor

% figure
%plot(time.*60./3.154e+7*365,heat1)
%xlabel('time [days]')
%ylabel('heat flux [W/m^2]')
%grid on
%grid minor

figure
plot(h_reference_signal,gload)
xlabel('altitude [km]')
ylabel('axial load factor [g]')
grid on
grid minor

figure
% plot(time.*60./3.154e+7,h)
% xlabel('time [years]')
plot(time_reference_signal.*60./3.154e+7*365,h_reference_signal)
xlabel('time [days]')
ylabel('altitude [km]')
grid on
grid minor

figure
% plot(time.*60./3.154e+7,V)
% xlabel('time [years]')
plot(time_reference_signal.*60./3.154e+7*365,h_reference_signal)
xlabel('time [days]')
ylabel('velocity [km/s]')
grid on
grid minor

figure
plot(gamma_reference_signal,h_reference_signal)
xlabel('flight path angle [deg]')
ylabel('altitude [km]')
grid on
grid minor

%Earth's surface circle generation
circ_ang=0:0.01:2.1*pi;
lcirc=length(circ_ang);
circ_r=ones(1,lcirc).*6371;

figure
polarplot(theta_reference_signal,h_reference_signal+6371,circ_ang,circ_r,theta_reference_signal(1),h_reference_signal(1)+6371,'*g',theta_reference_signal(end),h_reference_signal(end)+6371,'*r')
title('Trajectory shape evolution')
legend('s/c trajectory','Earth''s surface','Initial position','Landing position')
%Karman line crossing detection
kar_mask=fix(mean(find(h_reference_signal<100.1 & h_reference_signal>99.9)));
%isolation of the last orbit
up=theta_reference_signal(end)-2*pi+0.01;
down=theta_reference_signal(end)-2*pi-0.01;
orb_mask=fix(mean(find(theta_reference_signal<up & theta_reference_signal>down)));

figure
polarplot(theta_reference_signal(orb_mask:end),h_reference_signal(orb_mask:end)+6371,circ_ang,circ_r,theta_reference_signal(kar_mask),h_reference_signal(kar_mask)+6371,'or')
title('Re-entry trajectory')
legend('s/c trajectory','Earth''s surface','Karman line crossing')

figure
plot(Vx_reference_signal,h_reference_signal,Vz_reference_signal,h_reference_signal)
xlabel('velocity [km/s]')
ylabel('altitude [km]')
legend('Vx tangent','Vz radial')
grid on
grid minor