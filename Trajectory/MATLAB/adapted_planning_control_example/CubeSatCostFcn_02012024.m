function L = CubeSatCostFcn_02012024(stage, x, u, pvcost)
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
p  = pvcost(17);                         % Prediction Horizon



if stage == p + 1
    % Calculate cost for the terminal stage
    L = 0.5*(x - xf).'*Sf*(x - xf);
else
    % Calculate cost for intermediate stages
    L = 0.5*(x - xf).'*Q*(x - xf) + 0.5*u.'*R*u;
end
end
