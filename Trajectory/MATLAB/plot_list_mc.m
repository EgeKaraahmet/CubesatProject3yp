%% Plot
n_hx = size(List_h_reference_signal, 2);
%%%%%%%%%%
%% h vs. time
figure;
for i = 1:n_hx
    plot(List_time_reference_signal(:, i).*60./3.154e+7*365, List_h_reference_signal(:, i),'-');  
    hold on;  
end
xlabel('time [days]')
ylabel('altitude [km]')
grid on
grid minor
hold off;

%% V vs. h
figure;
for i = 1:n_hx
    plot(List_V_reference_signal(:, i), List_h_reference_signal(:, i),'-');  
    hold on;  
end
xlabel('velocity [km/s]')
ylabel('altitude [km]')
grid on
grid minor
hold off;

%% h vs. heat flux
figure;
for i = 1:n_hx
    plot(List_h_reference_signal(:, i), List_heat1(:, i),'-');  
    hold on;  
end
xlabel('altitude [km]')
ylabel('heat flux [W/m^2]')
grid on
grid minor
hold off

%%  heat against time 
figure;
for i = 1:n_hx
    plot(List_time_reference_signal(:, i).*60./3.154e+7*365, List_heat1(:, i),'-');  
    hold on;  
end

xlabel('time [days]')
ylabel('heat flux [W/m^2]')
grid on
grid minor
hold off
 %% h against t 
figure;
for i = 1:n_hx
    plot(List_time_reference_signal(:, i).*60./3.154e+7*365, List_h_reference_signal(:, i),'-');  
    hold on;  
end
xlabel('time [days]')
ylabel('altitude [km]')
grid on
grid minor
hold off

%% V against t
figure;
for i = 1:n_hx
    plot(List_time_reference_signal(:, i).*60./3.154e+7*365, List_V_reference_signal(:, i),'-');  
    hold on;  
end
xlabel('time [days]')
ylabel('velocity [km/s]')
grid on
grid minor
hold off


%% h vs. gload
figure;
for i = 1:n_hx
    plot(List_h_reference_signal(:, i), List_gload(:, i),'-');  
    hold on;  
end
xlabel('altitude [km]')
ylabel('axial load factor [g]')
grid on
grid minor

%% h vs gamma
figure;
for i = 1:n_hx
    plot(List_h_reference_signal(:, i), List_gamma_reference_signal(:, i),'-');  
    hold on;  
end
xlabel('flight path angle [deg]')
ylabel('altitude [km]')
grid on
grid minor

% %%
% %Earth's surface circle generation
% circ_ang=0:0.01:2.1*pi;
% lcirc=length(circ_ang);
% circ_r=ones(1,lcirc).*6371;
% 
% figure
% polarplot(theta_reference_signal,h_reference_signal+6371,circ_ang,circ_r,theta_reference_signal(1),h_reference_signal(1)+6371,'*g',theta_reference_signal(end),h_reference_signal(end)+6371,'*r')
% title('Trajectory shape evolution')
% legend('s/c trajectory','Earth''s surface','Initial position','Landing position')
% %Karman line crossing detection
% kar_mask=fix(mean(find(h_reference_signal<100.1 & h_reference_signal>99.9)));
% %isolation of the last orbit
% up=theta_reference_signal(end)-2*pi+0.01;
% down=theta_reference_signal(end)-2*pi-0.01;
% orb_mask=fix(mean(find(theta_reference_signal<up & theta_reference_signal>down)));
% % 
% figure
% polarplot(theta_reference_signal(orb_mask:end),h_reference_signal(orb_mask:end)+6371,circ_ang,circ_r,theta_reference_signal(kar_mask),h_reference_signal(kar_mask)+6371,'or')
% title('Re-entry trajectory')
% legend('s/c trajectory','Earth''s surface','Karman line crossing')
% 
% figure
% plot(Vx_reference_signal,h_reference_signal,Vz_reference_signal,h_reference_signal)
% xlabel('velocity [km/s]')
% ylabel('altitude [km]')
% legend('Vx tangent','Vz radial')
% grid on
% grid minor
