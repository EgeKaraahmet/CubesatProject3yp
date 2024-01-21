%% Land a Rocket Using Multistage Nonlinear MPC
% The Code is adapted from a MATLAB Example 

% In this planning and control problem:
% * position of the CubeSat is bounded in X (horizontal axis) from -150 to
% 150km. 
% * The goal position is at (0,0). 
% * The maximum thrust applied by each thruster can be preconfigured.
% * The rocket can have an arbitrary initial position and orientation. 

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

Ts = 1; 
%cpLander = size(h_reference_signal);
pLander = 1; % pLander = pLander(1)


%% Design Lander and Follow the Optimal Path
% Like generic nonlinear MPC, you can use multistage nonlinear MPC for
% reference tracking and disturbance rejection. In this example, you use it
% to track the optimal trajectory found by the planner. For a
% path-following problem, the lander does not require a long prediction
% horizon. Create the controller.

lander = nlmpcMultistage(pLander,5,1);
lander.Ts = Ts;

%%
% For the path-following controller, the lander has the same prediction
% model, thrust bounds, and minimum Y position.
lander.Model.StateFcn = 'CubeSatStateFcn_testing';
lander.Model.StateJacFcn = 'CubeSatStateJacobianFcn_testing';
lander.MV(1).Min = 0.2;
lander.MV(1).Max = 0.5;
lnader.States(1).Min = 0; 

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
% Since changes in control actions are represented by |MVRate|, you must
% enable the multistage nonlinear MPC controller to use |MVRate| values in
% its in calculations.
lander.UseMVRate = true;

%%
% Validate the controller design.
simdata = getSimulationData(lander);
validateFcns(lander,rand(5,1),rand(1,1),simdata);

%%
% Simulate the landing maneuver in a closed-loop control scenario by
% iteratively calling |nlmpcmove|. This simulation assumes that all states
% are measured. Stop the simulation when the rocket states are close
% enough to the target states.
x = x0;
u = u0;
k = 1;
references = reshape(info.Xopt',(pPlanner+1)*5,1); % Extract reference signal as column vector.
while true
    % Obtain new reference signals.
    simdata.StageParameter = RocketLanderReferenceSignal(k,references,pLander);
    % Compute the control action.
    [u,simdata,infoLander] = nlmpcmove(lander,x,u,simdata);
    % Update the animation plot.
    updatePlot(plotobj,(k-1)*Ts,x,u);    
    pause(0.1);
    % Simulate the plant to the next state using an ODE solver.
    [~,X] = ode45(@(t,x) RocketStateFcn(x,u),[0 Ts],x);
    x = X(end,:)';
    % Stop if rocket has landed.
    if max(abs(x-[0;10;0;0;0;0])) < 1e-2
        % Plot the the final rocket position.
        updatePlot(plotobj,k*Ts,x,zeros(2,1));    
        break
    end
    % Move to the next simulation step.
    k = k + 1;
end

%%
% Due to the shorter horizon and different cost terms, the landing
% trajectory slightly differs from the planned trajectory and it takes
% longer to land. This result is often what happens in such a two-tier
% control framework with planning and regulation.

%% Simulate in Simulink Using Multistage Nonlinear MPC Block
% You can implement the same closed-loop simulation in a Simulink model
% using the Multistage Nonlinear MPC block.
mdl = 'RocketLanderSimulation_testing';
open_system(mdl)

%%
% The States scope shows that the plant states are brought to the target
% states in a reasonable time.  
sim(mdl)
open_system([mdl '/States'])
open_system([mdl '/MVs'])

%%
% For real-time applications, you can generate code from the Multistage
% Nonlinear MPC block.

bdclose(mdl)
