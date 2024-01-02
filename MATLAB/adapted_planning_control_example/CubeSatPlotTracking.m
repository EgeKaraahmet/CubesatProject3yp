function CubeSatPlotTracking(Xopt, Ts, Psteps, Tsteps, Xcl, Ucl)
    tp = Ts * (0:Psteps);
    tt = Ts * (0:Tsteps);
    
    figure;
    states = {'h', 'X', 'V', 'theta', 'gamma'};
    
    for i = 1:5
        subplot(3, 2, i)
        plot(tt, Xcl(:, i), '+', tp, Xopt(:, i), '-')
        legend('actual', 'plan', 'location', 'northwest')
        title(states{i})
    end
    
    figure;
    for i = 1
        subplot(1, 1, i)
        stairs(tt(1:end-1), Ucl(:, i))
        title('Thrust u')
        axis([0 tt(end) -0.1 1.1])
        hold on
        stairs(tp(1:end-1), info.MVopt(1:end-1, i))
        legend('actual', 'plan')
        hold off
    end
    
    figure;
    hold on
    for ct = 1:size(Xopt, 1)
        lf = [cos(Xopt(ct, 4)) * 0.5, sin(Xopt(ct, 4)) * 0.5];
        rf = [cos(Xopt(ct, 4)) * 0.5, sin(Xopt(ct, 4)) * 0.5];
        lr = [-cos(Xopt(ct, 4)) * 0.5, -sin(Xopt(ct, 4)) * 0.5];
        rr = [-cos(Xopt(ct, 4)) * 0.5, -sin(Xopt(ct, 4)) * 0.5];
        patch([lf(1), rf(1), rr(1), lr(1)] + Xopt(ct, 2), ...
              [lf(2), rf(2), rr(2), lr(2)] + Xopt(ct, 1), 'y', 'FaceAlpha', 0.5, 'LineStyle', ':');
    end
    
    for ct = 1:size(Xcl, 1)
        lf = [cos(Xcl(ct, 4)) * 0.5, sin(Xcl(ct, 4)) * 0.5];
        rf = [cos(Xcl(ct, 4)) * 0.5, sin(Xcl(ct, 4)) * 0.5];
        lr = [-cos(Xcl(ct, 4)) * 0.5, -sin(Xcl(ct, 4)) * 0.5];
        rr = [-cos(Xcl(ct, 4)) * 0.5, -sin(Xcl(ct, 4)) * 0.5];
        if ct < size(Xcl, 1)
            patch([lf(1), rf(1), rr(1), lr(1)] + Xcl(ct, 2), ...
                  [lf(2), rf(2), rr(2), lr(2)] + Xcl(ct, 1), 'b', 'FaceAlpha', 0.5);
        else
            patch([lf(1), rf(1), rr(1), lr(1)] + Xcl(ct, 2), ...
                  [lf(2), rf(2), rr(2), lr(2)] + Xcl(ct, 1), 'r', 'FaceAlpha', 0.5);
        end
    end
    
    xlabel('X')
    ylabel('h')
    title('Compare with Optimal Trajectory')
    
    fprintf('Actual fuel consumption = %10.6f\n', sum(Ucl(1:end-1)) * Ts);
end
