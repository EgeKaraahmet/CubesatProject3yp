% %% Extract data
% load('A403.mat')
% h_plot_kmf_403 = xHistory_kmf(:,1);
% last_item = h_plot_kmf_403(end);  % Get the last item of the vector
% h_plot_kmf_403 = [h_plot_kmf_403; last_item];  % Append the last item to the vector
% X_plot_kmf_403 = xHistory_kmf(:,2);
% last_item = X_plot_kmf_403(end);  % Get the last item of the vector
% X_plot_kmf_403 = [X_plot_kmf_403; last_item];  % Append the last item to the vector
% 
% load('A402p5.mat')
% h_plot_kmf_402p5 = xHistory_kmf(:,1);
% X_plot_kmf_402p5 = xHistory_kmf(:,2);
% 
% load('A402.mat')
% h_plot_kmf_402 = xHistory_kmf(:,1);
% last_item = h_plot_kmf_402(end);  % Get the last item of the vector
% h_plot_kmf_402 = [h_plot_kmf_402; last_item];  % Append the last item to the vector
% X_plot_kmf_402 = xHistory_kmf(:,2);
% last_item = X_plot_kmf_402(end);  % Get the last item of the vector
% X_plot_kmf_402 = [X_plot_kmf_402; last_item];  % Append the last item to the vector
% 
% 
% load('A401p5.mat')
% h_plot_kmf_401p5 = xHistory_kmf(:,1);
% X_plot_kmf_401p5 = xHistory_kmf(:,2);
% 
% load('A401.mat')
% h_plot_kmf_401 = xHistory_kmf(:,1);
% X_plot_kmf_401 = xHistory_kmf(:,2);
% 
% load('A400p5.mat')
% h_plot_kmf_400p5 = xHistory_kmf(:,1);
% X_plot_kmf_400p5 = xHistory_kmf(:,2);
% 
% load('A400.mat')
% h_plot_kmf_400 = xHistory_kmf(:,1);
% X_plot_kmf_400 = xHistory_kmf(:,2);
% %%
% plot(abs(X_plot_kmf_400),h_plot_kmf_400, ...
%     abs(X_plot_kmf_400p5),h_plot_kmf_400p5, ...
%     abs(X_plot_kmf_401),h_plot_kmf_401, ...
%     abs(X_plot_kmf_401p5),h_plot_kmf_401p5, ...
%     abs(X_plot_kmf_402),h_plot_kmf_402, ...
%     abs(X_plot_kmf_402p5),h_plot_kmf_402p5, ...
%     abs(X_plot_kmf_403),h_plot_kmf_403, ...
%     abs(x_reference_atm),h_reference_atm,'b')
% 
% ylim([0, inf]);
% xlabel('Horizontal distance travelled (km)');
% ylabel('Altitude (km)');
% 
% 
% 
% %% 2cm 1mm
% x_ref_end = 1430210000; 
% x_kmf_end = x_ref_end * (1+0.1/20); 
% 
% h_alt = 6378; 
% 
% x_last_ref = x_ref_end - floor(x_ref_end/(2*pi*h_alt)); 
% x_last_kmf = x_kmf_end - floor(x_kmf_end/(2*pi*h_alt)); 
% 
% p_ref = sqrt(h_alt^2 + x_last_ref^2);
% p_kmf = sqrt(h_alt^2 + x_last_kmf^2);
% 
% abs(p_ref-p_kmf)

%% small area 
load('A122.mat')
h_plot_kmf_A122 = xHistory_kmf(:,1);
X_plot_kmf_A122 = xHistory_kmf(:,2);

load('A122p5.mat')
h_plot_kmf_A122p5 = xHistory_kmf(:,1);
X_plot_kmf_A122p5 = xHistory_kmf(:,2);

load('A123.mat')
h_plot_kmf_A123 = xHistory_kmf(:,1);
X_plot_kmf_A123 = xHistory_kmf(:,2);

load('A123p5.mat')
h_plot_kmf_A123p5 = xHistory_kmf(:,1);
X_plot_kmf_A123p5 = xHistory_kmf(:,2);

load('A124.mat')
h_plot_kmf_A124 = xHistory_kmf(:,1);
X_plot_kmf_A124 = xHistory_kmf(:,2);

load('A124p5.mat')
h_plot_kmf_A124p5 = xHistory_kmf(:,1);
X_plot_kmf_A124p5 = xHistory_kmf(:,2);

load('A125.mat')
h_plot_kmf_A125 = xHistory_kmf(:,1);
X_plot_kmf_A125 = xHistory_kmf(:,2);

load('A125p5.mat')
h_plot_kmf_A125p5 = xHistory_kmf(:,1);
X_plot_kmf_A125p5 = xHistory_kmf(:,2);

load('A126.mat')
h_plot_kmf_A126 = xHistory_kmf(:,1);
X_plot_kmf_A126 = xHistory_kmf(:,2);

load('A126p5.mat')
h_plot_kmf_A126p5 = xHistory_kmf(:,1);
X_plot_kmf_A126p5 = xHistory_kmf(:,2);

load('A127.mat')
h_plot_kmf_A127 = xHistory_kmf(:,1);
X_plot_kmf_A127 = xHistory_kmf(:,2);

load('A127p5.mat')
h_plot_kmf_A127p5 = xHistory_kmf(:,1);
X_plot_kmf_A127p5 = xHistory_kmf(:,2);

load('A128.mat')
h_plot_kmf_A128 = xHistory_kmf(:,1);
X_plot_kmf_A128 = xHistory_kmf(:,2);

%%
plot(abs(X_plot_kmf_A122)/1000,h_plot_kmf_A122/1000, ...
    abs(X_plot_kmf_A122p5)/1000,h_plot_kmf_A122p5/1000, ...
    abs(X_plot_kmf_A123)/1000,h_plot_kmf_A123/1000, ...
    abs(X_plot_kmf_A123p5)/1000,h_plot_kmf_A123p5/1000, ...
    abs(X_plot_kmf_A124)/1000,h_plot_kmf_A124/1000, ...
    abs(X_plot_kmf_A124p5)/1000,h_plot_kmf_A124p5/1000, ...
    abs(X_plot_kmf_A125)/1000,h_plot_kmf_A125/1000, ...
    abs(X_plot_kmf_A125p5)/1000,h_plot_kmf_A125p5/1000, ...
    abs(X_plot_kmf_A126)/1000,h_plot_kmf_A126/1000, ...
    abs(X_plot_kmf_A126p5)/1000,h_plot_kmf_A126p5/1000, ...
    abs(X_plot_kmf_A127)/1000,h_plot_kmf_A127/1000, ...
    abs(X_plot_kmf_A127p5)/1000,h_plot_kmf_A127p5/1000, ...
    abs(X_plot_kmf_A128)/1000,h_plot_kmf_A128/1000,...
    abs(x_reference_atm)/1000,h_reference_atm/1000,'b');
 

ylim([0, inf]);
xlabel('Horizontal distance travelled (km)');
title('Perturbation test - varying fin areas: altitude vs. horizontal distance travelled plot')
ylabel('Altitude (km)');
legend('Reference altitude vs. distance ', 'Trajectory Tracking (with EKF)')

% Add legend
legend('A=22cm^2','A=22.5cm^2','A=23cm^2', 'A=23.5cm^2', 'A=24cm^2', 'A=24.5cm^2', 'A=25cm^2','A=25.5cm^2','A=26cm^2','A=26.5cm^2','A=27cm^2','A=27.5cm^2','A=28cm^2',  'Reference Altitude', 'Location', 'southwest');
% Calculation 
h_alt = 6378; 
%%
X_last_ref = abs(x_reference_atm(end))/1000; 
%X_last_ref = X_last_ref - (2*pi*h_alt)*floor(X_last_ref/(2*pi*h_alt));
p_ref = sqrt(h_alt^2 + X_last_ref^2);


X_last_123 = abs(X_plot_kmf_A123(end))/1000; 
%X_last_123 = X_last_123 - (2*pi*h_alt)*floor(X_last_123/(2*pi*h_alt));
p_123 = sqrt(h_alt^2 + X_last_123^2); 
pe_123 = abs(p_123 - p_ref);

X_last_123p5 = abs(X_plot_kmf_A123p5(end))/1000; 
%X_last_123p5 = X_last_123p5 - (2*pi*h_alt)*floor(X_last_123p5/(2*pi*h_alt));
p_123p5 = sqrt(h_alt^2 + X_last_123p5^2);
pe_123p5 = abs(p_123p5 - p_ref);

X_last_124 = abs(X_plot_kmf_A124(end))/1000; 
%X_last_124 = X_last_124 - (2*pi*h_alt)*floor(X_last_124/(2*pi*h_alt));
p_124 = sqrt(h_alt^2 + X_last_124^2);
pe_124 = abs(p_124 - p_ref);

X_last_124p5 = abs(X_plot_kmf_A124p5(end))/1000; 
%X_last_124p5 = X_last_124p5 - (2*pi*h_alt)*floor(X_last_124p5/(2*pi*h_alt));
p_124p5 = sqrt(h_alt^2 + X_last_124p5^2);
pe_124p5 = abs(p_124p5 - p_ref);

X_last_125 = abs(X_plot_kmf_A125(end))/1000; 
%X_last_125 = X_last_125 - (2*pi*h_alt)*floor(X_last_125/(2*pi*h_alt));
p_125 = sqrt(h_alt^2 + X_last_125^2); 
pe_125 = abs(p_125 - p_ref);

X_last_125p5 = abs(X_plot_kmf_A125p5(end))/1000; 
%X_last_125p5 = X_last_125p5 - (2*pi*h_alt)*floor(X_last_125p5/(2*pi*h_alt));
p_125p5 = sqrt(h_alt^2 + X_last_125p5^2);
pe_125p5 = abs(p_125p5 - p_ref);
% 
X_last_126 = abs(X_plot_kmf_A126(end))/1000; 
%X_last_126 = X_last_126 - (2*pi*h_alt)*floor(X_last_126/(2*pi*h_alt));
p_126 = sqrt(h_alt^2 + X_last_126^2);
pe_126 = abs(p_126 - p_ref);

X_last_126p5 = abs(X_plot_kmf_A126p5(end))/1000; 
%X_last_126p5 = X_last_126p5 - (2*pi*h_alt)*floor(X_last_126p5/(2*pi*h_alt));
p_126p5 = sqrt(h_alt^2 + X_last_126p5^2);
pe_126p5 = abs(p_126p5 - p_ref);

X_last_127 = abs(X_plot_kmf_A127(end))/1000; 
%X_last_127 = X_last_127 - (2*pi*h_alt)*floor(X_last_127/(2*pi*h_alt));
p_127 = sqrt(h_alt^2 + X_last_127^2);
pe_127 = abs(p_127 - p_ref);

X_last_127p5 = abs(X_plot_kmf_A127p5(end))/1000; 
%X_last_127p5 = X_last_127p5 - (2*pi*h_alt)*floor(X_last_127p5/(2*pi*h_alt));
p_127p5 = sqrt(h_alt^2 + X_last_127p5^2);
pe_127p5 = abs(p_127p5 - p_ref);

X_last_128 = abs(X_plot_kmf_A128(end))/1000; 
%X_last_128 = X_last_128 - (2*pi*h_alt)*floor(X_last_128/(2*pi*h_alt));
p_128 = sqrt(h_alt^2 + X_last_128^2);
pe_128 = abs(p_128 - p_ref);

% Define the range of areas in cm^2
areas = [23, 23.5, 24, 24.5, 25, 25.5, 26, 26.5, 27, 27.5, 28];

% Define the pe values
pe = [pe_123, pe_123p5, pe_124, pe_124p5, pe_125,pe_125p5, pe_126, pe_126p5,pe_127,pe_127p5, pe_128 ];%, pe_125p5, pe_126];
pe = pe * 100; 
% Create a table
T = table(areas', pe', 'VariableNames', {'Area (cm^2)', 'Position Error (km)'});

% Display the table
disp(T);

%%
plot(areas,pe,'*-')
xlabel('Perturbation Area (cm^2)');
title('Perturbation test - Position Error')
ylabel('Position Error (km)')

