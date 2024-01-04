clear 
close all 
clc


%% Initial conditions 
% Input 
A=100*10^(-4);          %m^2        % spacecraft cross - sectional area

% constant 
m = 6;                  %kg         % spacecraft mass  m = 3 for ISS
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
h0=kepOut.h(end);            %m        %height
x0=kepOut.x(end);            %m        %travelled distance
V0=kepOut.V(end);            %m/s      %orbital speed
theta0=kepOut.theta(end);    %rad      %true anomaly
gamma0=kepOut.gamma(end);    %rad      %flight path angle


X0 = x0; 
%Re-entry simulation, atmospheric phase, more accuracy is needed.
%Smaller integration step required (e.g. 0.1 s)
%Simulation stops when altitude reaches 0 km (default)
%%%%%%%%%%%%%%%%%%
stop_h=0;   %m          %final desired altitude
%%%%%%%%%%%%%%%%%%

atmOut = sim('reference_signal_generator','Solver',solv_atm,'FixedStep',step_atm,'StartTime','kepOut.time(end)');

h_end=atmOut.h(end);         %m        %altitude
x_end=atmOut.x(end);         %m        %travelled distance
V_end=atmOut.V(end);            %m/s      %orbital speed
theta_end=atmOut.theta(end);    %rad      %true anomaly
gamma_end=atmOut.gamma(end);    %rad      %flight path angle


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





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Define Objective Cost Function Initial Values 
% x: (1) h: altitude 
%    (2) X: distance travelled
%    (3) V: velocity 
%    (4) theta: true anomaly
%    (5) gamma: flight path angle
% Stage Parameters

Sf = [20; 20; 0; 0; 0];     % Terminal Weight Matrix.
Q  = [30; 30; 0; 0; 0];     % State Weight Matrix.
R  = 10;                   % Control Weighting Matrix.


p  = 40;                    % Prediction Horizon.
Ts = 10000;                  % Sampling time 

xf = [h_end; x_end; V_end; theta_end; gamma_end];       % Terminal State.

% Combine stage parameters into a column array.
pvcost = [xf; Sf; Q; R; p];

% Define initial conditions
x_op_0 = [h0;X0;V0;theta0;gamma0]; 
x_op_0 = [200*10^3;0;7788;0;0]; 



%% Design Nonlinear MPC Controller
% Create a nonlinear multistage MPC object with 5 states and 1 inputs. 
nx  = 5;
nmv = 1;
msobj = nlmpcMultistage(p,nx,nmv);

msobj.Ts = Ts;

%% Specify the prediction model state using the CubeSat dynamics functions.
msobj.Model.StateFcn = "CubeSatStateFcn_25122023";


%% Specify this cost function using a named function
% Setting cost function for each control interval.
for k = 1:p+1
    msobj.Stages(k).CostFcn = "CubeSatCostFcn_02012024";
    msobj.Stages(k).ParameterLength = length(pvcost);
end

%% Specify hard bounds. 
msobj.MV.Min = 100 * 10^(-4);
msobj.MV.Max =  150 * 10^(-4);



%% Initialize data structure
% The state and stage functions require state and stage parameters. 
% Use getSimulationData to initialize data structure. 
simdata = getSimulationData(msobj);
% simdata.StateFcnParameter = pvstate;
simdata.StageParameter = repmat(pvcost, p+1, 1);

%% Setting the Optimization Solver: conjugate gradient method
msobj.Optimization.Solver = "cgmres";
% Adjust the Stabilization Parameter 
% based on the prediction model sample time.
msobj.Optimization.SolverOptions.StabilizationParameter = 1/msobj.Ts; % 1/mosbj.Ts

% Set the solver parameters.
msobj.Optimization.SolverOptions.MaxIterations = 10;
msobj.Optimization.SolverOptions.Restart = 3;
msobj.Optimization.SolverOptions.BarrierParameter = 1e-3;
msobj.Optimization.SolverOptions.TerminationTolerance = 1e-6;


%% Simulation duration in seconds.
Duration = Ts * p;

x0 = x_op_0; 
u0 = 100 * 10^(-4); 


% Store states and control for plotting purposes.
xHistory1 = x0.';
uHistory1 = u0.';

% Initialize control.
uk = uHistory1(1,:);

% Initialize the accumulated elapsed time 
% for computing the optimal control action calculation.
timerVal1 = 0;

% Simulation loop
for k = 1:(Duration/Ts)
    % Compute optimal control action using nlmpcmove.
    xk = xHistory1(k,:).';
    
    % Call the nlmpcmove function.
    tic
    [uk, simdata] = nlmpcmove(msobj, xk, uk, simdata);
    
    % Accumulate the elapsed time.
    timerVal1 = timerVal1 + toc;

    % Simulate CubeSat trajectory for the next control interval.
    ODEFUN = @(t,xk) CubeSatStateFcn_25122023(xk,uk);
    [TOUT, XOUT] = ode45(ODEFUN, [0 Ts], xHistory1(k,:));

    % Log states and control.
    xHistory1(k+1,:) = XOUT(end,:);
    uHistory1(k+1,:) = uk;
   

    % Check for clicked Cancel button.
    %if getappdata(hbar,"canceling")
    %    delete(hbar)
    %    break
    % end
end


%%
% Assuming h_plot is xHistory1(:,1)
h_plot = xHistory1(:,1);
X_plot = xHistory1(:,2);
T_plot = TOUT;
u_plot = uHistory1; 

% Create a logical index to exclude negative values
positive_indices = h_plot >= 0;

% Filter the vectors using the logical index
h_plot = h_plot(positive_indices);
X_plot = X_plot(positive_indices);
T_plot = T_plot(positive_indices);
u_plot = u_plot(positive_indices);



%%
figure
plot(X_plot,h_plot,'r')
figure
plot(x_reference_signal,h_reference_signal,'b')

