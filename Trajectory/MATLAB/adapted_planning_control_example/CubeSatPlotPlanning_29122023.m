function CubeSatPlotPlanning_29122023(Info,Ts)
    % CubeSatPlotPlanning displays the optimal trajectory of the flying
    % cubesat 
    
    Xopt = Info.Xopt;
    MVopt = Info.MVopt;
    fprintf('Minimum cost function value = %10.6f\n', Info.Cost * Ts)
    
    % Examine solution
    t = Info.Topt;
    figure;
    states = {'h', 'X', 'V', 'theta', 'gamma'};
    for i = 1:size(Xopt, 2)
        subplot(3, 2, i)
        plot(t, Xopt(:, i), 'o-')
        title(states{i})
    end
    
    % Replace the last row u(k+p) with 0
    MVopt(end, :) = 0;
    
    figure;
    stairs(t, MVopt, 'o-')
    title('Control Input u');
    
    % Plot optimal trajectory with modified X values
    figure;
    X_modified = Xopt(:, 2);
    % X_modified(Xopt(:, 2) > 0) = -abs(Xopt(Xopt(:, 2) > 0, 2));
    plot(X_modified, Xopt(:, 1), 'o-') % plot h against X
    xlabel('X')
    ylabel('h')
    title('Optimal Trajectory')
end
