% ReadMe.m 
%%
%
% % State x
% % x: (1) h: altitude 
% %    (2) X: distance travelled
% %    (3) V: velocity 
% %    (4) theta: true anomaly
% %    (5) gamma: flight path angle
% 
% % input: 
% % u: fin areas 
% 
% This folder contains the following files:
% 1. main_03012024.m:  
%   the main script that contains a reference state generator, 
%   followed by a nonlinear controller design. Please run this file.
% 
%       - initial state: 
%           x_op_0 = [h0;X0;V0;theta0;gamma0]; generated from the reference state generator 
%       - final desired state: 
%           xf = [0; x_end; V_end; theta_end; gamma_end];generated from the reference state generator 
% 
%
%
% 2. CubeSatStateFcn_25122023.m: 
%       the nonlinear system, taking 
%                   input states: x =[h;X;V;theta; gamma] and input: u, 
%                   output state first derivative: dxdt. 
% 
% 
% 3. CubeSatCostFcn_03122024: cost function containing the following terms: 
%
%% This cost function is adapted from the MATLAB example: Control Robot Manipulator Using C/GMRES Solver
%   -> you may refer to the example by using the command: 
openExample('mpc/ControlRobotManipulatorUsingCGMRESSolverExample')
%% Explanation of the cost functions 
% 
% *    •    Stage cost: The stage cost represents the cost associated with each stage (the time step) 
%                       of the control problem. It quantifies the immediate cost incurred at each step. 
%                       In this example, the stage cost is defined using a standard 
%                       quadratic cost function that penalizes the
%                       deviation of the state vector x from the desired state x_f at each time point within the horizon,k, 
%                       and the control effort u, required to achieve the desired state. 
%
%                       I have also included another cost penalty to
%                       penalise the maximum heat flux in the stage cost
%  
% *    •    Terminal cost: This term represents a cost associated with the final state of the system. 
%                           The terminal penalty penalizes the deviation of the state at the end of the horizon x_p+1, 
%                           from the desired state x_f, where p+1 is the final stage, 
%                           which corresponds to the end of the prediction horizon,t+T .  
% 
% Here, the matrices Q and R determine the weight of the state deviation from the desired state, and the control effort, respectively. 
% The matrix Sf determines the weight for the terminal penalty. 
% The overall performance index for the receding horizon control problem is
% a sum of the stage cost and terminal cost. 
% 
% 4.   reference_signal_generator.slx
%      This is a simulink code used for reference state generator 
% 
% 
% 
