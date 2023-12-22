function J = RocketLanderCostFcn_testing(stage,x,u,dmv,h_reference_signal)
% Rocket lander cost function.

% Copyright 2020 The MathWorks, Inc.



if stage == 1
    J = dmv'*[0.1 0;0 0.1]*dmv;
elseif stage == 11
    J = (x-h_reference_signal)'*(x-h_reference_signal);
else
    J = (x-h_reference_signal)'*(x-h_reference_signal) + dmv'*[stage*0.1 0;0 stage*0.1]*dmv;
end