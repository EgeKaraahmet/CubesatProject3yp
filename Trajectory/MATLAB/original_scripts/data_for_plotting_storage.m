
%%
V_my_model = V; 
x_my_model = x; 
h_my_model = h;
heat1_my_model = heat1; 
gload_my_model = gload; 
time_my_model = time; 
gamma_my_model = gamma; 
theta_my_model = theta; 

ndays_my_model=ndays;
gload_max_my_model = gload_max;
heat1_max_my_model = heat1_max; 

%%
%% Plotting

% Figure 1
figure
plot(V, h, 'b', V_my_model, h_my_model, 'r')
xlabel('Velocity [km/s]')
ylabel('Altitude [km]')
grid on
grid minor
legend('De Cecio''s Work', 'Group Design')
title('Trajectory Dynamics System - Velocity vs Altitude')

% Figure 2
figure
plot(-x, h, 'b', -x_my_model, h_my_model, 'r')
xlabel('Horizontal Distance tavelled [km]')
ylabel('Altitude [km]')
grid on
grid minor
legend('De Cecio''s Work', 'Group Design')
title('Trajectory Dynamics System - Horizontal Distance vs Altitude')

%
% Filter out the parts with h >= 150 km
% Filter out the parts with h >= 150 km
filtered_indices = h < 100;
filtered_x = -x(filtered_indices);
filtered_h = h(filtered_indices);

% Use the filtered indices for h to filter h_my_model
filtered_indices = h_my_model < 100;
filtered_h_my_model = h_my_model(filtered_indices);
filtered_x_my_model = x_my_model(filtered_indices);

% Plot the filtered data
figure;
plot(filtered_x, filtered_h, 'b', -filtered_x_my_model, filtered_h_my_model, 'r');
xlabel('Horizontal Distance traveled [km]');
ylabel('Altitude [km]');
grid on;
grid minor;
legend('De Cecio''s Work', 'Group Design');
title('Trajectory Dynamics System - Horizontal Distance vs Altitude (Altitude < 150 km)');



% Figure 3
figure
plot(h, heat1, 'b', h_my_model, heat1_my_model, 'r')
xlabel('Altitude [km]')
ylabel('Heat Flux [W/m^2]')
grid on
grid minor
legend('De Cecio''s Work', 'Group Design')
title('Trajectory Dynamics System - Altitude vs Heat Flux')

% Find and label the maximum point on the curves
[maxHeat1, idxHeat1] = max(heat1);
text(h(idxHeat1), maxHeat1, sprintf(' Max: %.2f W/m^2', maxHeat1), 'Color', 'b', 'FontSize', 10, 'HorizontalAlignment', 'right');

[maxHeat1_my_model, idxHeat1_my_model] = max(heat1_my_model);
text(h_my_model(idxHeat1_my_model), maxHeat1_my_model, sprintf(' Max: %.2f W/m^2', maxHeat1_my_model), 'Color', 'r', 'FontSize', 10, 'HorizontalAlignment', 'left');

% Figure 4
figure
plot(h, gload, 'b', h_my_model, gload_my_model, 'r')
xlabel('Altitude [km]')
ylabel('Axial Load Factor [g]')
grid on
grid minor
legend('De Cecio''s Work', 'Group Design')
title('Trajectory Dynamics System - Altitude vs Axial Load Factor')

% Find and label the maximum point on the curves
[maxGload, idxGload] = max(gload);
text(h(idxGload), maxGload, sprintf(' Max: %.2f g', maxGload), 'Color', 'b', 'FontSize', 10, 'HorizontalAlignment', 'right');

[maxGload_my_model, idxGload_my_model] = max(gload_my_model);
text(h_my_model(idxGload_my_model), maxGload_my_model, sprintf(' Max: %.2f g', maxGload_my_model), 'Color', 'r', 'FontSize', 10, 'HorizontalAlignment', 'left');

% Figure 5
figure
plot(time .* 60 ./ 3.154e+7 * 365, h, 'b', time_my_model .* 60 ./ 3.154e+7 * 365, h_my_model, 'r')
xlabel('Time [days]')
ylabel('Altitude [km]')
grid on
grid minor
legend('De Cecio''s Work', 'Group Design')
title('Trajectory Dynamics System - Time vs Altitude')


%% error cal
abs(heat1_max_my_model-heat1_max)/heat1_max*100

abs(gload_max_my_model-gload_max)/gload_max*100

abs(ndays_my_model-ndays)/ndays*100

abs(x_my_model(end) - x(end)) / abs(x(end)) * 100


alti = (200+6371); 
%%
no_of_orbit_my_model = floor(x_my_model(end)/(2*pi*6371));
x_distance_my_model = x_my_model(end)-no_of_orbit_my_model*(2*pi*6371);
distance_my_model = sqrt(alti^2 + (x_distance_my_model)^2) - 6371

no_of_orbit = floor(x(end)/(2*pi*6371));
x_distance = x(end)-no_of_orbit*(2*pi*6371);
distance = sqrt(alti^2 + (x_distance)^2) - 6371
