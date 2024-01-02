function CubeSatPlotPlanning_25122023(Info,Ts)
% CubeSatPlotPlanning displays the optimal trajectory of the flying
% robot.

% Copyright 2018-2021 The MathWorks, Inc.

Xopt = Info.Xopt;
MVopt = Info.MVopt;
fprintf('Minimum cost function value = %10.6f\n',Info.Cost*Ts)

%% Examine solution
t = Info.Topt;
figure;
states = {'h','X','V','theta','gamma'};
for i = 1:size(Xopt,2)
    subplot(3,2,i)
    plot(t,Xopt(:,i),'o-')
    title(states{i})
end
figure;
MVopt(end,:) = 0; % replace the last row u(k+p) with 0
stairs(t,MVopt,'o-')
title('Control Input u');

figure;
plot(Xopt(:,2),Xopt(:,1),'o-') % plot h against X
xlabel('X')
ylabel('h')
title('Optimal Trajectory')