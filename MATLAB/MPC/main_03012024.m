clear
close all
%% 
% Initial A = 135 cm^2,    Radius = 150 km 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Reference state generator 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initial conditions 
% Input 
A=110*10^(-4);          %m^2        % spacecraft cross - sectional area

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

kepOut = sim('reference_signal_generator_MPC','Solver',solv_kep,'FixedStep',step_kep);

%   updating initial conditions, 2nd run using ouput from 1st

%space environment
ndays1=kepOut.time(end)/60/60/24;  %days since intial time of simulation
jdate=jdate+ndays1;                %new julian date at beginning of 2nd run
F107_avg=90.85;         %SFU: updating
F107_day=86.0;          %SFU: updating
Kp=1;

space_coes = [jdate;F107_avg;F107_day;Kp];

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

atmOut = sim('reference_signal_generator_MPC','Solver',solv_atm,'FixedStep',step_atm,'StartTime','kepOut.time(end)');

h_end=atmOut.h(end);         %m        %altitude
x_end=atmOut.x(end);         %m        %travelled distance
V_end=atmOut.V(end);            %m/s      %orbital speed
theta_end=atmOut.theta(end);    %rad      %true anomaly
gamma_end=atmOut.gamma(end);    %rad      %flight path angle


%% Results analysis
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

testing = atmOut.time; 
x_state_reference = [atmOut.h, atmOut.x, atmOut.V, atmOut.theta, atmOut.gamma];




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% nonlinear MPC Controller 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Define Objective Cost Function Initial Values 
% State x
% x: (1) h: altitude 
%    (2) X: distance travelled
%    (3) V: velocity 
%    (4) theta: true anomaly
%    (5) gamma: flight path angle

% input: 
% u: fin areas 

% Stage Parameters

Sf = [100; 100; 100; 100; 100];     % Terminal Weight Matrix.
Q  = [30; 30; 30; 30; 30];          % State Weight Matrix.
R  = 10;                            % Control Weighting Matrix.
Q_heat = 100;                       % Heat Flux Weight Matrix


% MPC controller setup
%%%% duration = 3851
p  = 4000/(100);                    % Prediction Horizon.
Ts = 100;                  % Sampling time 

xf = [h_end; x_end; V_end; theta_end; gamma_end];       % Desired terminal State.


% xf = [target_h; target_X; target_V; target_theta; target_gamma];
% Combine stage parameters into a column array.
pvcost = [xf; Sf; Q; R; Q_heat; p;space_coes];



% Define initial conditions
x_op_0 = [h0;X0;V0;theta0;gamma0]; 
% x_op_0 = [200*10^3;0;7788;0;0]; 




%% Design Nonlinear MPC Controller
% Create a nonlinear multistage MPC object with 5 states and 1 inputs. 
nx  = 5;
nmv = 1;
msobj = nlmpcMultistage(p,nx,nmv);

msobj.Ts = Ts;

%% Specify the prediction model state using the CubeSat dynamics functions.
msobj.Model.StateFcn = "CubeSatStateFcn_25122023";
msobj.Model.StateJacFcn = "CubeSatStateJacobianFcn_25122023"; 


%% Specify this cost function using a named function
%% Setting cost function for each control interval.
for k = 1:p+1
    msobj.Stages(k).CostFcn = "CubeSatCostFcn_03012024";
     msobj.Stages(k).CostJacFcn = "CubeSatCostJacobianFcn_03012024";
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

%% Setting the Optimization Solver: C/GMRES
%%% continuation/generalized minimum residual (C/GMRES) solver
msobj.Optimization.Solver = "cgmres";
% Adjust the Stabilization Parameter 
% based on the prediction model sample time.
msobj.Optimization.SolverOptions.StabilizationParameter = 1/msobj.Ts; % 1/mosbj.Ts

% Set the solver parameters.
msobj.Optimization.SolverOptions.MaxIterations = 100;
msobj.Optimization.SolverOptions.Restart = 3;
msobj.Optimization.SolverOptions.BarrierParameter = 1e-5;
msobj.Optimization.SolverOptions.TerminationTolerance = 1e-6;




%% Simulation duration in seconds.
Duration = Ts * p;

x0 = x_op_0; 
u0 = A; 


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

    % Compute and store the cost for the current iteration
    L = CubeSatCostFcn_03012024(p,xHistory1(k,:), uk,pvcost)
end



%% 
% Assuming h_plot is xHistory1(:,1)
h_plot = xHistory1(:,1);
X_plot = xHistory1(:,2);
V_plot = xHistory1(:,3);
theta_plot = xHistory1(:,4);
gamma_plot = xHistory1(:,5);
T_plot = TOUT;
u_plot = uHistory1; 

% % % Create a logical index to exclude negative values
positive_indices = h_plot >= 0;
% % 
% % % Filter the vectors using the logical index
h_plot = h_plot(positive_indices);
X_plot = X_plot(positive_indices);
V_plot = V_plot(positive_indices);
theta_plot = theta_plot(positive_indices);
gamma_plot = gamma_plot(positive_indices);
% T_plot = T_plot(positive_indices);
u_plot = u_plot(positive_indices);

x_plot_out = [h_plot X_plot V_plot theta_plot gamma_plot];

h_reference_atm = [atmOut.h];
x_reference_atm = [atmOut.x];

[p_op,~] = size(h_plot);
%%
% Plotting using zero-order hold
figure
subplot(2,2,1)
n = 1:length(u_plot);
stairs(n, u_plot, 'b-', 'LineWidth', 2);

% Adding labels and title
ylabel('u');
title('Zero-Order Hold Plot');

%%%
% I doubt that the controller is not working properly. This is because my
% u-plot is almost a step increase to 125 cm^2, even when I change my
% sampling time. 


subplot(2,2,2)
plot(X_plot, h_plot, 'r', x_reference_atm, h_reference_atm, 'b')
title('Nonlinear MPC output')

% Adding legend
legend('Nonlinear MPC output (red)', 'Reference altitude vs. distance (blue)')


A_target = u_plot(end);   %m^2        % spacecraft cross - sectional area
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Feedback Control and Kalman Filter 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Feedback Control for Path Following
% After the optimal trajectory is found, a feedback controller is required
% to move the robot along the path. In theory, you can apply the optimal MV
% profile directly to the thrusters to implement feed-forward control.
% However, in practice, a feedback controller is needed to reject
% disturbances and compensate for modeling errors.
%
% You can use different feedback control techniques for tracking. In this
% example, you use a generic nonlinear MPC controller to move the robot to
% the final location. In this path tracking problem, you track references
% for all 5 states (the number of outputs equals the number of states).
ny = 5;
msobj_tracking = nlmpc(nx,ny,nmv);

%%
% Use the same state function and its Jacobian function.
msobj_tracking.Model.StateFcn = msobj.Model.StateFcn;

%%
% For tracking control applications, reduce the computational effort by
% specifying shorter prediction horizon (no need to look far into the
% future) and control horizon (for example, free moves are allocated at the
% first few prediction steps).
msobj_tracking.Ts = Ts;
msobj_tracking.PredictionHorizon = 10;
msobj_tracking.ControlHorizon = 4;

%%
% The default cost function in nonlinear MPC is a standard quadratic cost
% function suitable for reference tracking and disturbance rejection. For
% tracking, tracking error has higher priority (larger penalty weights on
% outputs) than control efforts (smaller penalty weights on MV rates).
msobj_tracking.Weights.ManipulatedVariablesRate = 10; 
msobj_tracking.Weights.OutputVariables = 100*ones(1,ny);

%%
% Set the same bounds for the area inputs.
msobj_tracking.MV.Min = 100 * 10^(-4);
msobj_tracking.MV.Max = 150 * 10^(-4);


%%
validateFcns(msobj_tracking,x0,u0);

%% Nonlinear State Estimation
% In this example, only the position states (h,X) are
% measured. The velocity states are unmeasured and must be estimated. Use
% an extended Kalman filter (EKF) from Control System Toolbox(TM) for
% nonlinear state estimation.
%
% Because an EKF requires a discrete-time model, you use the trapezoidal
% rule to transition from x(k) to x(k+1), which requires the solution of
% |nx| nonlinear algebraic equations. For more information, open
% |FlyingRobotStateFcnDiscreteTime.m|.
DStateFcn = @(xk,uk,Ts) CubeSatStateFcnDiscreteTime_25122023(xk,uk,Ts);

%%
% Measurement can help the EKF correct its state estimation
DMeasFcn = @(xk) xk(1:2);

%%
% Create the EKF, and indicate that the measurements have noise.
noise_coe = 1; 
EKF = extendedKalmanFilter(DStateFcn,DMeasFcn,x0);
EKF.MeasurementNoise = noise_coe;


%% Closed-Loop Simulation of Tracking Control
% Simulate the system for 32 steps with correct initial conditions.
Tsteps = p_op;
xHistory_kmf = x0';
uHistory_kmf = [];

lastMV = A_target;
%%
% The reference signals are the optimal state trajectories computed at the
% planning stage. When passing these trajectories to the nonlinear MPC
% controller, the current and future trajectory is available for
% previewing.
% Xopt_kmf = x_plot_out;
%[p_kmf,~] = size(Xopt_kmf);
% Xref_kmf = [Xopt_kmf(1:p_kmf,:);repmat(Xopt_kmf(end,:),Tsteps-p_kmf,1)];
Xref_kmf = x_plot_out; 
% Xref_kmf = x_state_reference; 

%%
% Use |nlmpcmove| and |nlmpcmoveopt| command for closed-loop simulation.
hbar = waitbar(0,'Simulation Progress');
options = nlmpcmoveopt;

for k = 1:Tsteps

    % Obtain plant output measurements with sensor noise.
    yk = xHistory_kmf(k,1:2)' + randn*noise_coe;

    % Correct state estimation based on the measurements.
    xk = correct(EKF, yk);

    % Compute the control moves with reference previewing.
    [uk,options] = nlmpcmove(msobj_tracking,xk,lastMV,Xref_kmf(k:min(k+9,Tsteps),:),[],options);

    % Predict the state for the next step.
    predict(EKF,uk,Ts);

    % Store the control move and update the last MV for the next step.
    uHistory_kmf(k,:) = uk'; %#ok<*SAGROW>
    lastMV = uk;

    % Update the real plant states for the next step by solving the
    % continuous-time ODEs based on current states xk and input uk.
    ODEFUN = @(t,xk) CubeSatStateFcn_25122023(xk,uk);
    [TOUT,YOUT] = ode45(ODEFUN,[0 Ts], xHistory_kmf(k,:)');

    % Store the state values.
    xHistory_kmf(k+1,:) = YOUT(end,:);

    % Update the status bar.
    waitbar(k/Tsteps, hbar);
end
close(hbar)

%%
% Assuming h_plot is xHistory1(:,1)
% Assuming you have the following vectors already defined
h_plot_kmf = xHistory_kmf(:,1);
X_plot_kmf = xHistory_kmf(:,2);
V_plot_kmf = xHistory_kmf(:,3);
theta_plot_kmf = xHistory_kmf(:,4);
gamma_plot_kmf = xHistory_kmf(:,5);
u_plot_kmf = uHistory_kmf;

% Create a logical index to exclude negative values
positive_indices_kmf = h_plot_kmf >= 0;

% Filter the vectors using the logical index
h_plot_kmf_non_negative = h_plot_kmf(positive_indices_kmf);
X_plot_kmf_non_negative = X_plot_kmf(positive_indices_kmf);
V_plot_kmf_non_negative = V_plot_kmf(positive_indices_kmf);
theta_plot_kmf_non_negative = theta_plot_kmf(positive_indices_kmf);
gamma_plot_kmf_non_negative = gamma_plot_kmf(positive_indices_kmf);
% u_plot_kmf_non_negative = u_plot_kmf(positive_indices_kmf);

% Find the greatest negative element for h_plot_kmf
h_kmf_negative_indices = h_plot_kmf < 0;
h_kmf_greatestNegativeElement = max(h_plot_kmf(h_kmf_negative_indices));
index_greatestNegativeElement = find(h_plot_kmf == h_kmf_greatestNegativeElement);

% Store the greatest negative element in all other vectors
h_plot_kmf = [h_plot_kmf_non_negative; h_kmf_greatestNegativeElement];
X_plot_kmf = [X_plot_kmf_non_negative; X_plot_kmf(index_greatestNegativeElement)];
V_plot_kmf = [V_plot_kmf_non_negative; V_plot_kmf(index_greatestNegativeElement)];
theta_plot_kmf = [theta_plot_kmf_non_negative; theta_plot_kmf(index_greatestNegativeElement)];
gamma_plot_kmf = [gamma_plot_kmf_non_negative; gamma_plot_kmf(index_greatestNegativeElement)];
% u_plot_kmf = [u_plot_kmf_non_negative; u_plot_kmf(h_kmf_negative_indices)];



%%
subplot(2,2,3)
plot(X_plot, h_plot, 'r',X_plot_kmf,h_plot_kmf,'*g')
title('OP, KMF')

subplot(2,2,4)
plot(X_plot, h_plot, 'r', x_reference_atm, h_reference_atm, 'b',X_plot_kmf,h_plot_kmf,'*g')
title('RS, OP, KMF')

% Adding legend
legend('Nonlinear MPC output (red)', 'Reference altitude vs. distance (blue)','kmf (green)')

