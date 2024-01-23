%% Land a Rocket Using Multistage Nonlinear MPC
%%
% States:
%
% x: (1) h: altitude 
%    (2) X: distance travelled
%    (3) V: velocity 
%    (4) theta: true anomaly
%    (5) gamma: flight path angle
%
% Input: 
% u: fin areas
%  

%% Initialisation 
x0 = [200*10^3;0;7788;0;0];  
u0 = 0.5;
Ts = 0.2;
pPlanner = 50;

%%
% Create a multistage nonlinear MPC controller for the specified prediction
% horizon and sample time.
planner = nlmpcMultistage(pPlanner,5,1);
planner.Ts = Ts;

%%
% Specify prediction model and its analytical Jacobian.
planner.Model.StateFcn = 'CubeSatStateFcn';
planner.Model.StateJacFcn = 'CubeSatStateJacobianFcn';

%% 
% Specify hard bounds on the two thrusters. You can adjust the maximum
% thrust and observe its impact on the landing strategy chosen by the
% planner. 
planner.MV(1).Min = 0;
planner.MV(1).Max = 0.5; 

%% 
% To avoid crashing, specify a hard lower bound on the vertical Y position.
planner.States(1).Min = 0;

%% 

for ct=1:pPlanner
    planner.Stages(ct).CostFcn = 'RocketPlannerCostFcn_testing';
    planner.Stages(ct).CostJacFcn = 'RocketPlannerCostGradientFcn_testing';
end

%% 
% To ensure a successful landing at the target, specify terminal state for
% the final stage. 
planner.Model.TerminalState = [0;0;0;0;0];

%% 
% In this example, set the maximum number of iterations to a large value
% to accommodate the large search space and the nonideal default initial
% guess.
planner.Optimization.SolverOptions.MaxIterations = 1000;

%% 
% After creating you nonlinear MPC controller, check whether there is any
% problem with your state, cost, and constraint functions, as well as their
% analytical Jacobian functions. To do so, call |validateFcns| functions
% with random initial plant states and inputs.
validateFcns(planner,rand(5,1),rand(1,1));

%%
% Compute the optimal landing path using |nlmpcmove|, which can typically
% take a few seconds, depending on the initial rocket position.
fprintf('Rocker landing planner running...\n');
tic;
[~,~,info] = nlmpcmove(planner,x0,u0);
t=toc;
fprintf('Calculation Time = %s\n',num2str(t));
fprintf('Objective cost = %s',num2str(info.Cost));
fprintf('ExitFlag = %s',num2str(info.Iterations));
fprintf('Iterations = %s\n',num2str(info.Iterations));

%%
% Extract the optimal trajectory from the |info| structure and plot the
% result.
figure
subplot(2,1,1)
plot(info.Xopt(:,1),info.Xopt(:,2),'*')
title('Optimal XY Trajectory')
% subplot(2,1,2)
% plot(info.Topt,info.MVopt(:,1),info.Topt,info.MVopt(:,2))
% title('Optimal MV Trajectory')



%% Design Lander and Follow the Optimal Path
% Like generic nonlinear MPC, you can use multistage nonlinear MPC for
% reference tracking and disturbance rejection. In this example, you use it
% to track the optimal trajectory found by the planner. For a
% path-following problem, the lander does not require a long prediction
% horizon. Create the controller.
pLander = 10;
lander = nlmpcMultistage(pLander,5,1);
lander.Ts = Ts;

%%
% For the path-following controller, the lander has the same prediction
% model, thrust bounds, and minimum Y position.
lander.Model.StateFcn = 'CubeSatStateFcn';
lander.Model.StateJacFcn = 'CubeSatStateJacobianFcn';
lander.MV(1).Min = 0;
lander.MV(1).Max = 0.5;
planner.States(1).Min = 0;
%%
% The cost function for the lander is different from that of the planner.
% The lander uses quadratic cost terms to achieve both tight reference
% tracking (by penalizing the tracking error) and smooth control actions
% (by penalizing large changes in the control actions). This lander cost
% function is implemented in the |RocketLanderCostFcn| function. The
% corresponding manually derived cost gradient function is implemented in
% |RocketLanderCostGradientFcn|.
%
% At run time, you provide the six state trajectory references to the
% lander as stage parameters. Therefore, specify the number of parameters
% for each stage.

for ct=1:pLander+1
    lander.Stages(ct).CostFcn = 'RocketLanderCostFcn_testing';
    lander.Stages(ct).CostJacFcn = 'RocketLanderCostGradientFcn_testing';
    lander.Stages(ct).ParameterLength = 5;
end

%%
