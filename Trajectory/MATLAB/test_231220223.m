% 定义系统动力学模型
sys = 'controller_dynamics'; % 替换为你的系统模型名称

% 定义性能指标
% performance = @(A) CostFcn_23122023(A);

% 设置约束条件
lb = 0.1;  % Lower bound for input A
ub = 0.5;  % Upper bound for input A

% 选择优化算法
algorithm = 'ga'; % 替换为你选择的优化算法

% 创建优化选项
options = optimoptions(algorithm, 'ConstraintTolerance', 1e-6, 'MaxGenerations', 500);

% 定义优化问题
optimal_trajectory = ga(sys, 1, [], [], [], [], lb, ub, [], options);

% 分析结果
optimal_trajectory


function J = CostFcn_23122023(A)
% cost function.
J = sum(A);
end