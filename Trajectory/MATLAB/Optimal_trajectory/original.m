%=========================================================================%
% PARAM_OPT.m %
%
% This is the parameter optimization program used in Chapter 6. %
% The first half of the code (prior to the opt. loop) can be used to
% run a single re-entry simulaion.
% %
% Derrick Tetzman
% April 2010 %=========================================================================%
clear 
close all 
clc
format long

%Planet-specific constants: Earth entry. 
rho_0 = 1.225; %kg/m^3.
H = 6.93E3; %Scale height in m. 
g_0 = 9.81; %Grav acc m/s^2.
R = 6378000; %Radius of Earth in m.
C = 7.28E-4; %Heat flux constant as adjusted previously.
a = 308.4; %Speed of sound in m/s (constant with isothermal assumption).

%Initial conditions 
h_0 = 121.92E3; %m 
x_0 = 0;
V_0 = 11137; %m/s 
q_0 = 0;

%Spacecraft constants: Apollo 4.
B = 372; %mass/CD*A ballistic coefficient in kg/m^2 
Rn = 4.661; % nose radius of vehicle in m
m = 5357; % mass of command module in kg

%% Read in Apollo 4 data and wind data from text files
% fileID = fopen('Ap4LD.txt'); 
% uAP4 = fscanf(fileID,'%f'); 
% fclose(fileID);
% usim = uAP4;%%%%%%%%%%%%%%%%% This sets the initial control input to % the approximated values for Apollo 4.
% Make sure to set n = 19 if using this.
% load('wind_profile.mat')

%% Control conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
deg = 6.93; % Desired initial flight path angle.
% Apollo 4 initial flight path angle = 6.93 deg from horiz. % Taken Directly from the Apollo heating data for AS-501.
gamma_0 = deg*pi/180; %Convert to Radians

x_target = 3797; % %Desired Final Range in km, AP4 final range 3797 km 

beta = 10^-11; %Step size--> 10^-11 used for ch. 6 tests.

n = 3; %Number of normalized time intervals. Apollo 4 = 19
Nmax = 1000; %Max number of iterations.
Imin = 1; %Max value of I for convergence.

delta_u(:,1) = -0.2*ones(n,1); %Initial step
t_f = 589.5; %Sensible guess for initial final time,
%based on previous simulation studies.
wind = 0; %Use if adding a constant wind along with block in Simulink
% ... units: m/s. Positive wind = Positive X direction.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%
%% Constraints
% % %
% Max heating and loading allowed, based on the max heating
% ... and loading from Apollo 4. Max loading experienced on Apollo 4 ... was roughly L_max = 7.25 g and max heat flux was q_max = 425 ... BTU/ft^2-s or 482.65 W/cm^2 or 4.83*10^6 W/m^2.
% ... See Sep 16th tech memo.
q_limit = 795.96*10^6; %max heat flux allowed in w/m^2 
L_limit = 12; %max loading allowed in g's.
h_limit = h_0; %Max altitude allowed in m.

%% First step in the algorithm
% Initialize variables
e = zeros(n,1); % 1xn zero vector for "interval" block, used because
                % embedded m-code doesn't let you reset the vector % within the block for whatever reason.

Id(:,1) = ones(n,1); %Placeholder for the first element of Id. 
cc = zeros(1,n); %initial constraint vector value.
vio = 0; %used to identify if a constraint has be violated.

% BUILD INITIAL GUESS VECTOR FOR CONTROL INPUT
for jj=1:1:n
    usim(jj,1) = 0.2; %Constant initial Guess.
    % if jj >= n-1
    % usim(jj,1) = 0.3 - jj/10; %Ascending / Descending Guess
end
%usim = [0.35 -0.3 0.2 0.2 0.2 0.2]'; % Custom initial guess. 
u(:,1) = usim;

% INTERVAL DEFINITION LOOP
for j=1:1:n
    tau(j) = j/n*t_f;
end % time interval definition

sim('ENTRYSIM'); %This runs the Simulink model using above values. 
maxX(1) = max(x);
maxT(1) = max(time);
maxH(1) = max(h);
maxQ(1) = max(q_dot); 
maxL(1) = max(L);

% CONSTRAINT CHECK LOOP
for row=1:1:length(c) 
    for col=1:1:n
        if row > 1
            if c(row,col) > c(row-1,col) 
                cc(col) = c(row,col);
                vio = 1;
            end
        end

    end 
end

% COMPUTE OVERALL COST FUNCTION
I(1) = (maxX - x_target)^2;

if vio == 1 %Stop if initial guess violates constraints
    disp('Initial guess violates constraints. Try another guess.') 
    disp(['Constraint vector = ',num2str(cc)])
    fprintf('1 = heating \n2 = loading \n3 = altitude') 
    fprintf('4 = heating + loading \n5 = heating + altitude')
    fprintf('6 = loading + altitude \n7 = All constraints \n')
elseif I(1) < Imin
    disp('Initial guess minimizes cost function. Done.')
else
% INTERVAL DEFINITION UPDATE 
    t_f = maxT(1); %Capture final time
    for j=1:1:n
        tau(j) = j/n*t_f;
    end

    % First Parameter update
    u(:,2) = usim + delta_u(:,1); %initial step is a sensible constant 
    last = delta_u(:,1);

    %% Optimization Loop-------------------------------------------------- 
    tic %Measure performance time of opt. loop
    
    for i=2:1:Nmax
    
        usim = u(:,i); %for input into multisim 
        disp(num2str(usim')) 
        sim('ENTRYSIM');
        maxH(i) = max(h);
        maxQ(i) = max(q_dot); 
        maxL(i) = max(L); 
        maxX(i) = max(x); 
        maxT(i) = max(time);
    
        % INTERVAL DEFINITION UPDATE 
        t_f = maxT(i); %Capture final time
        for j=1:1:n
            tau(j) = j/n*t_f;
        end
        % COMPUTE OVERALL COST FUNCTION
        I(i) = (maxX(i) - x_target)^2; 
        disp(num2str(I(i)))
        if I(i) < Imin
            disp('I < Imin. Close enough. Terminating loop.')
            break 
        end
        cc = zeros(1,n);
        vio = 0;
    
        % CONSTRAINT CHECK LOOP 
        for row=1:1:length(c)
            for col=1:1:n
            if row == 1
                if c(row,col) > 0 
                   cc(col) = c(row,col); 
                   vio = 1;
                end
            elseif row > 1
                if c(row,col) > c(row-1,col) 
                cc(col) = c(row,col);
                vio = 1;
                end 
            end
            end 
        end

        if vio == 1 
            disp(num2str(cc))
        end
    
        disp(' ')
        if i == Nmax %value reached after Nmax iterations 
            disp('Max number of iterations reached. i = Nmax.') 
            break % Further update unnecessary.
        end
    
    
        % DERIVATIVE / PARAMETER UPDATE LOOP ------
        t_f = maxT(i-1); % Use previous final time 
        for j=1:1:n
            tau(j) = j/n*t_f;
        end


        for k=1:1:n
            if (cc(k) ~= 0) %Constraint Control
            % ADD CONSTRAINT LOGIC HERE %
            % EXAMPLE BRUTE FORCE METHOD I USED:
                if (cc(k)==1) || (cc(k)==2) || (cc(k)==4) 
                    u(k,i+1) = u(k,i) + 0.1;
                    if k ~= 1
                        u(k-1,i+1) = u(k-1,i+1) + 0.0001 ;
                    end
                elseif cc(k)==3
                    u(k,i+1) = u(k,i) - 0.1 ;
                    if k ~= 1
                        u(k-1,i+1) = u(k-1,i+1) - 0.001 ; 
                    end
        
                else
                    u(k,i+1) = u(k,i-1);
                end
             Id(k,i) = 10^-6; %Placeholder for the derivative array
    
            elseif (u(k,i) == u(k,i-1)) 
                u(k,i+1) = u(k,i) - 10^-3; 
                disp('NO UPDATE!') 
                disp(['interval = ',num2str(k)]) 
                beep
            %Alerts the user if algorithm gets stuck.
            %This typically happens when the spacecraft skips out.
            else
                usim = u(:,i-1);
                usim(k) = u(k,i);
                sim('ENTRYSIM')
                maxsX = max(x);
                Id(k,i) = ((maxsX - x_target)^2-I(i-1)) / (usim(k)-u(k,i-1));
                if Id(k,i) == 0 
                    Id(k,i) = 10^-4;
                end
                % PARAMETER UPDATE
                delta_u(k,i) = -beta * Id(k,i); 
                u(k,i+1) = u(k,i) + delta_u(k,i); 
                last = delta_u(:,i);
            end 
        end
            t_f = maxT(i); %Switch back for the next loop iteration. 
            for j=1:1:n
                tau(j) = j/n*t_f;
            end
            % Note: it might seem strange to have these additional loops 
            % ...for tau just for the derivative loop, but it saves memory 
            % ...by eliminating the need to store yet ANOTHER large array.
    end


%% Display Results MODIFIED FOR ONLY ONE RUN THROUGH!
disp('************* RESULTS **************') 
disp(' ')
disp(['For n = ',num2str(n)])
disp(' ')
if vio == 1
    disp('Constraints violated.')
    disp(['Constraint vector = ', num2str(cc)])
    fprintf('1 = heating \n2 = loading \n3 = altitude')
    fprintf('4 = heating + loading \n5 = heating + altitude')
    fprintf('6 = loading + altitude \n7 = All constraints \n')
else
    disp('No constraint violations.') 
end

toc

disp(' ')
disp(['# Iterations: ',num2str(i)])
disp(' ')
disp('The optimum values of L/D are:') 
disp(num2str(u(:,i)))
disp(' ')
disp('The final cost function value is:') 
disp(['I = ',num2str(I(i))])
disp(' ')
disp('With a corresponding final range of:') 
disp(['x_final = ',num2str(maxX(i))])

%% Read in Approximate AP4 Trajectory XY reconstruction
fileIDX = fopen('AP4apxrecX.txt'); 
fileIDY = fopen('AP4apxrecY.txt'); 
Xrec = fscanf(fileIDX,'%f');
Yrec = fscanf(fileIDY,'%f'); fclose(fileIDX);
fclose(fileIDY);

%% Plot
% Plot Final Trajectory
figure5 = figure('Color',[1 1 1]);
hold on
plot(x,h,'b','DisplayName','Simulated Trajectory'); 
plot(Xrec,Yrec,'k','DisplayName','Approx. Reconstruction'); 
grid on; 
xlabel('Range (km)'); 
ylabel('Altitude (km)'); 
title('Final Trajectory');
hold off
end % For the initial constraint violation check.

beep % Alerts user that program has completed. %% "Interval" embedded M-function for reference.
... This chooses the control input given the current time.
% function [L_D,c] = interval(f,n,usim,tau,t,q_limit,L_limit,q_dot,L,h,h_limit)
% %#eml %
% %MAKE SURE L IS ABS VALUE OF LOADING
% L_D = usim(1); %
% d = f; %
%
% for k = 1:1:n %
% if t <= tau(k)
% L_D = usim(k); %
% if (q_dot >= q_limit) && (h <= h_limit/1000) && (L < L_limit)
% d(k) = 1;
% elseif (L >= L_limit) && (h <= h_limit/1000) && (q_dot < q_limit)
% d(k) = 2;
% elseif (h > h_limit/1000) && (q_dot < q_limit) && (L < L_limit)
% d(k) = 3;
% elseif (q_dot >= q_limit) && (L >= L_limit) && (h <= h_limit/1000)
% d(k) = 4;
% elseif (q_dot >= q_limit) && (h > h_limit/1000) && (L < L_limit)
% d(k) = 5;
% elseif (L >= L_limit) && (h > h_limit/1000) && (q_dot < q_limit)
% d(k) = 6;
% elseif (L >= L_limit) && (h > h_limit/1000) && (q_dot >= q_limit)
% d(k) = 7;
% end
% break
% end %
% end
% c = d;
% end
%% "hwindv" embedded M-function for reference.
% function wind = hwindv(h,windEast)
% %#eml %
% if h > 75
% wind = 0;
% else
% i = round(h);
% wind = windEast(i);
% end %
% end