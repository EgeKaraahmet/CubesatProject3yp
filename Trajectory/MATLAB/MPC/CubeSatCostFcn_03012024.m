function L = CubeSatCostFcn_03012024(stage, x, u, pvcost)
% twolinkCostFcn
%   L = twolinkCostFcn(stage, x, u, pvcost)
%   
%   This function computes the objective cost for the two-link robot
%   manipulator example.
%
%   Inputs:
%       - stage     Current stage of the Multistage MPC problem
%       - x         Current state vector
%       - u         Control action vector at the previous interval
%       - pvcost    Vector containing various stage parameters
%
%   Output:
%       - L:        Cost value calculated by the cost function.

% Copyright 2023 The MathWorks, Inc.


% Extract Terminal Desired State
xf = pvcost(1:5);

% Extract and Initialize Wiegthing Matrices
Sf = diag(pvcost(6:10));                 % Terminal Weight Matrix
Q  = diag(pvcost(11:15));                % State Weight Matrix
R  = diag(pvcost(15+(1:length(u))));     % Control Weight Matrix
Q_heat = pvcost(17);
p  = pvcost(end);                         % Prediction Horizon

h = x(1);
V = x(3);
%q_ref = 200 * 10^3; 
rcurv = 0.1; 

% heat flux
SH = 8397.5; 
rho = 1.225 * exp(- h/(SH));    % my Python model, using SH = 8.43 and rho0=1.221 as constant

q_max=1.83e-4*V.^3.*sqrt(rho/rcurv);
q_max_reference = 200 * 10^3; 



if stage == p + 1
    % Define the distance constraint parameters
    % gamma_mpc = 1; % Adjust this value based on your requirements
    % V = diag([100, 100, 100, 100, 100]); 

    % Calculate cost for the terminal stage with additional distance constraint
    % distance_term = 0.5 * (x - xf).' * V * (x - xf);
    % terminal_cost = 0.5 * (x - xf).' * Sf * (x - xf);
    % L = terminal_cost + gamma_mpc * distance_term;
    L = 0.5 * (x - xf).' * Sf * (x - xf);
else
    % Calculate cost for intermediate stages
    % L = 0.5*(x - xf).'*Q*(x - xf) + 0.5*u.'*R*u + (q_max - q_max_reference) * Q_heat *(q_max-q_max_reference);
    % L = 0.5*(x - xf).'*Q*(x - xf) + (q_max - q_ref) * Q_heat *(q_max-q_ref);
    L = 0.5*(x - xf).'*Q*(x - xf) + 0.5*u.'*R*u; 
end
end
