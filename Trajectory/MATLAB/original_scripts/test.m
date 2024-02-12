% Define a range of heights for testing
test_heights = linspace(0, 210, 100); % Generate 100 test heights from 0 to 210 km

% Initialize arrays to store rho and SH values
rho_values = zeros(size(test_heights));
SH_values = zeros(size(test_heights));

% Call the rho_ege function for each test height
for i = 1:length(test_heights)
    [rho_values(i), SH_values(i)] = rho_ege(test_heights(i));
end

% Plot the results
figure;
subplot(2, 1, 1);
plot(test_heights, rho_values, 'b-', 'LineWidth', 2);
xlabel('Height (km)');
ylabel('\rho (kg/m^3)');
title('Density vs. Height');
grid on;

subplot(2, 1, 2);
plot(test_heights, SH_values, 'r-', 'LineWidth', 2);
xlabel('Height (km)');
ylabel('SH');
title('SH vs. Height');
grid on;
