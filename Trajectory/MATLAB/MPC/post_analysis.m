%% Plots 
close all
%% plots 
%%
% Plotting using zero-order hold
figure
subplot(1,2,1)
n = 1:length(u_plot_op);
stairs(n, u_plot_op, 'b-', 'LineWidth', 2);


% Adding labels and title
ylabel('u');
title('Zero-Order Hold Plot');

subplot(1,2,2)
n = 1:length(u_plot_kmf);
stairs(n, u_plot_kmf, 'b-', 'LineWidth', 2);


% Adding labels and title
ylabel('u');
title('Zero-Order Hold Plot');
%%
figure

n = 1:length(u_plot_op);
stairs(n, u_plot_op - 100*10^(-4), 'b-', 'LineWidth', 2);

% Adding labels and title
xlabel('Stages');
ylabel('Aerobrake Areas (cm^2)');
title('Optimal Control Law (U*)');






%%
figure

subplot(1,3,1)
plot(abs(X_plot_op)/1000, h_plot_op/1000, 'r', abs(x_reference_atm)/1000, h_reference_atm/1000, 'b')
title('RS , OP')

% Adding legend
legend('Trajectory Tracking (no EKF)', 'Reference altitude vs. distance ')

xlabel('Horizontal distance travelled (km)');
ylabel('Altitude (km)');

subplot(1,3,2)
plot(abs(X_plot_op)/1000, h_plot_op/1000, 'r',abs(X_plot_kmf)/1000,h_plot_kmf/1000,'*g')
title('OP, KMF')

xlabel('Horizontal distance travelled (km)');
ylabel('Altitude (km)');

subplot(1,3,3)
plot(abs(X_plot_op)/1000, h_plot_op/1000, 'r', abs(x_reference_atm)/1000, h_reference_atm/1000, 'b',abs(X_plot_kmf)/1000,h_plot_kmf/1000,'*g')
title('RS, OP, KMF')

% Adding legend
legend('Trajectory Tracking (no EKF)', 'Reference altitude vs. distance ','Trajectory Tracking (with EKF)')

xlabel('Horizontal distance travelled (km)');
ylabel('Altitude (km)');


%%
% Plot 1: Nonlinear MPC output
figure
plot(abs(X_plot_op)/1000, h_plot_op/1000, 'r', abs(x_reference_atm)/1000, h_reference_atm/1000, 'b')
title('Altitude vs. horizontal distance travelled plot')
xlabel('Horizontal distance travelled (km)');
ylabel('Altitude (km)');
legend('Trajectory Tracking (no EKF)', 'Reference altitude vs. distance ')

% Plot 2: OP, KMF
figure
plot(abs(X_plot_op)/1000, h_plot_op/1000, 'r', abs(X_plot_kmf)/1000, h_plot_kmf/1000, '*g')
title('Altitude vs. horizontal distance travelled plot')
ylim([0, inf]);
xlabel('Horizontal distance travelled (km)');
ylabel('Altitude (km)');
legend('Trajectory Tracking (no EKF)', 'Trajectory Tracking (with EKF)')

% Plot 3: RS, KMF
figure
plot(abs(x_reference_atm)/1000, h_reference_atm/1000, 'b', abs(X_plot_kmf)/1000, h_plot_kmf/1000, '-*g')
title('Altitude vs. horizontal distance travelled plot')
ylim([0, inf]);
xlabel('Horizontal distance travelled (km)');
ylabel('Altitude (km)');
legend('Reference altitude vs. distance ', 'Trajectory Tracking (with EKF)')

% Plot 4: RS, OP, KMF
figure
plot(abs(X_plot_op)/1000, h_plot_op/1000, 'r', abs(x_reference_atm)/1000, h_reference_atm/1000, 'b', abs(X_plot_kmf)/1000, h_plot_kmf/1000, '-*g')
title('Altitude vs. horizontal distance travelled plot')
ylim([0, inf]);
xlabel('Horizontal distance travelled (km)');
ylabel('Altitude (km)');
legend('Trajectory Tracking (no EKF)', 'Reference altitude vs. distance ', 'Trajectory Tracking (with EKF)')




%% Calculations 
altitude = 200 + 6378; 

x_ref_landing = abs(x_reference_atm(end))/1000;
x_kmf_landing = abs(X_plot_kmf(end))/1000;
position_ref = sqrt(altitude^2+x_ref_landing^2);
position_kmf = sqrt(altitude^2+x_kmf_landing^2);

position_kmf- position_ref