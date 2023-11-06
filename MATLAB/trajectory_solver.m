%% Constants 
deg    = math.pi/180;   % convert degrees to radians
g0     = 9.81;          %  sea-level acceleration of gravity (m/s)
Re     = 6378e3;        % radius of the earth (m)
rho0   = 1.225;         % sea level density of atmosphere (kg/m^3)

R = 287.4;              % specific gas constant for air
K = 1.4;                % ratio of specific heat of a gas at a constant pressure to heat at a constant volume (1.4 for air)
T0 = 288.15;            % K standard temperature at sea level

% Cubesat values
Rn     = 0.25;          % nose radius (m)
m      = 5;             % Mass at orbit (kg)

% aerodynamic constant
S      = 0.785;         % frontal area, m^2
k      = 2;           % aerodynamic constant, 0.03




classdef Planet
    properties
        radius
        grav_k
        name
    end

    methods
        function obj = Planet(radius, grav_k, name)
            obj.radius = radius;
            obj.grav_k = grav_k;
            obj.name = name;
        end

        function g = gravity(obj, r)
            g = -(obj.grav_k / (r ^ 2.0));
        end

        function alt = altitude(obj, pos)
            alt = norm(pos) - obj.radius;
        end

        function cart = cartesian(obj, lat, alt)
            r = alt + obj.radius;
            cart = [r * cos(lat), r * sin(lat)];
        end

        function pol = polar(obj, x, y)
            alt = sqrt(x ^ 2 + y ^ 2) - obj.radius;
            pol = [atan2(x, y), alt];
        end

        function dens = density(~, alt)
            dens = 0;  % Replace with the actual density function if needed
        end
    end
end

classdef Earth < Planet
    methods
        function obj = Earth()
            obj = obj@Planet(6371e3, 3.986004418e14, 'Earth');
        end

        function dens = density(obj, alt)
            dens = 1.221 * exp(-alt / 8.43e3);
        end
    end
end

function [x, y, vx, vy, ax, ay, t] = sim_run(sim, planet, craft)
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

function do_plot(xlabel, x, ylabel, y, label, title, fname)
    figure;
    plot(x, y, 'DisplayName', label);
    title(title);
    xlabel(xlabel);
    ylabel(ylabel);
    legend;
    grid on;
    tightfig;
    % Save the figure if needed: saveas(gcf, fname, 'png');
end

bodies.Earth = Earth;

% Load parameters from a JSON file (replace with your parameter file path)
params = jsondecode(fileread('parameter_file.json'));

craft = params.craft;
planet = bodies.(params.planet);
sim = params.sim;

% Run the simulation
[x, y, vx, vy, ax, ay, t] = sim_run(sim, planet, craft);

title = sprintf('%s (\\beta=%.0f kg/m^2) @ %s', craft.name, craft.ballistic_coef, planet.name);
label = sprintf('fpa=%.2f, v_EI=%.0f m/s', sim.fpa, sim.velocity);

% Convert to useful coordinates
[lat, alt] = planet.polar(x, y);
downrange = lat * planet.radius;
vn = norm([vx, vy]);

% Compute the velocity magnitude
v = sqrt(vx .^ 2 + vy .^ 2);

% Get the axial load ((ax,ay) projected onto (vx,vy))
aa = abs((ax .* vx + ay .* vy) / v);

% Time and distance to go
tti = max(t) - t;
dtg = max(downrange) - downrange;

f1 = do_plot('downrange (km)', downrange / 1e3, 'altitude (km)', alt / 1e3, label, title, sprintf('%s-traj.png', craft.name));
f2 = do_plot('axial loads (g)', aa / 9.81, 'altitude (km)', alt / 1e3, label, title, sprintf('%s-load_alt.png', craft.name));
f3 = do_plot('time since EI (s)', t, 'axial loads (g)', aa / 9.81, label, title, sprintf('%s-load_time.png', craft.name));
f4 = do_plot('distance to splashdown (km)', dtg / 1e3, 'time to parachute deploy (s)', tti / 60, label, title, sprintf('%s-dtg.png', craft.name));
f5 = do_plot('velocity (km/s)', v / 1e3, 'altitude (km)', alt / 1e3, label, title, sprintf('%s-vel.png', craft.name));
