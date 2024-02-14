%% Extract data
load('A403.mat')
h_plot_kmf_403 = xHistory_kmf(:,1);
last_item = h_plot_kmf_403(end);  % Get the last item of the vector
h_plot_kmf_403 = [h_plot_kmf_403; last_item];  % Append the last item to the vector
X_plot_kmf_403 = xHistory_kmf(:,2);
last_item = X_plot_kmf_403(end);  % Get the last item of the vector
X_plot_kmf_403 = [X_plot_kmf_403; last_item];  % Append the last item to the vector

load('A402p5.mat')
h_plot_kmf_402p5 = xHistory_kmf(:,1);
X_plot_kmf_402p5 = xHistory_kmf(:,2);

load('A402.mat')
h_plot_kmf_402 = xHistory_kmf(:,1);
last_item = h_plot_kmf_402(end);  % Get the last item of the vector
h_plot_kmf_402 = [h_plot_kmf_402; last_item];  % Append the last item to the vector
X_plot_kmf_402 = xHistory_kmf(:,2);
last_item = X_plot_kmf_402(end);  % Get the last item of the vector
X_plot_kmf_402 = [X_plot_kmf_402; last_item];  % Append the last item to the vector


load('A401p5.mat')
h_plot_kmf_401p5 = xHistory_kmf(:,1);
X_plot_kmf_401p5 = xHistory_kmf(:,2);

load('A401.mat')
h_plot_kmf_401 = xHistory_kmf(:,1);
X_plot_kmf_401 = xHistory_kmf(:,2);

load('A400p5.mat')
h_plot_kmf_400p5 = xHistory_kmf(:,1);
X_plot_kmf_400p5 = xHistory_kmf(:,2);

load('A400.mat')
h_plot_kmf_400 = xHistory_kmf(:,1);
X_plot_kmf_400 = xHistory_kmf(:,2);

plot(abs(X_plot_kmf_400),h_plot_kmf_400, ...
    abs(X_plot_kmf_400p5),h_plot_kmf_400p5, ...
    abs(X_plot_kmf_401),h_plot_kmf_401, ...
    abs(X_plot_kmf_401p5),h_plot_kmf_401p5, ...
    abs(X_plot_kmf_402),h_plot_kmf_402, ...
    abs(X_plot_kmf_402p5),h_plot_kmf_402p5, ...
    abs(X_plot_kmf_403),h_plot_kmf_403, ...
    abs(x_reference_atm),h_reference_atm,'b')

ylim([0, inf]);
xlabel('Horizontal distance travelled (km)');
ylabel('Altitude (km)');
legend('Trajectory Tracking (no EKF)', 'Reference altitude vs. distance ', 'Trajectory Tracking (with EKF)')
