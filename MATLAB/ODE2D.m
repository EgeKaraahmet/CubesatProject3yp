function [x, y, vx, vy, ax, ay, t] = ODE2D(sim, planet, craft)
    max_it = sim.max_it;
    dt = sim.delta_t;
    fpa = deg2rad(sim.fpa);

    p = [0, sim.entry_interface + planet.radius];
    x = zeros(1, max_it);
    y = zeros(1, max_it);

    v = sim.velocity * [cos(fpa), sin(fpa)];
    vx = zeros(1, max_it);
    vy = zeros(1, max_it);

    a = [0, 0];
    ax = zeros(1, max_it);
    ay = zeros(1, max_it);

    t = 0:dt:(max_it - 1) * dt;

    beta = craft.ballistic_coef;
    ld = craft.lift_drag;

    k = 1;
    for i = 1:max_it
        p = p + v * dt;
        x(k) = p(1);
        y(k) = p(2);

        r = norm(p);
        rho = planet.density(planet.altitude(p));
        v_mag = norm(v);
        normal = [v(2), v(1)];

        aero_accel = 0.5 * rho * v_mag * (ld * normal / beta - v / beta);
        gravity_accel = planet.gravity(r) * (p / r);

        a = aero_accel + gravity_accel;
        ax(k) = a(1);
        ay(k) = a(2);

        v = v + a * dt;
        vx(k) = v(1);
        vy(k) = v(2);

        k = k + 1;

        if planet.altitude(p) <= sim.stop_alt
            fprintf('done in %d iterations\n', k);
            break;
        end
    end

    x = x(1:k - 1);
    y = y(1:k - 1);
    vx = vx(1:k - 1);
    vy = vy(1:k - 1);
    ax = ax(1:k - 1);
    ay = ay(1:k - 1);
    t = t(1:k - 1);
end

