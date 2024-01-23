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

%% Initial conditions 
% Input 
A=(100) * (10)^(-2*2);   %m^2        % spacecraft cross - sectional area

% constant 
m = 6;                    %kg         % spacecraft mass  m = 3 for ISS
Cd=2.2;                             % drag coefficient
BC = m / Cd / A; 



%constants
mi=398600.44;           %km^3/s^2   % Earth G*M
Re=6371.0088;           %km         % Earth mean radius
g0=9.80665;             %m/s^2      % Gravitational acceleration (sea level)
%spacecraft parameters
E=0;                                % aerodynamic efficiency
CL=Cd*E;                            % lift coefficient
rcurv=0.1;              %m          % TPS nose curvature radius

%space environment      
date=[2023,01,01,12,00,00];         % intial time of simulation [y,m,d,h,m,s]
jdate=juliandate(date);
F107_avg=90.85;         %SFU        % F10.7 average of 3x27 days before the date under consideration
F107_day=86.0;          %SFU        % F10.7 average of day before the date under consideration


Kp=1;                                 % Kp three-hourly planetary geomagneticindex

alt_0 = 200;              %km         % initial altitude    %% 400 km
r_a=alt_0+Re;             %km         % radius at apoapsis
r_p=alt_0+Re;             %km         % radius at periapsis
a=(r_a+r_p)/2;            %km         % semimajor axis
e=(r_a-r_p)/(r_a+r_p);                % eccentricity


%re-entry path initial values (at t=0 s)
x0=0;                             % m          % travelled distance
gamma0=0;                         % rad        % flight path angle
theta0=0;                         % rad        % true anomaly 
r0=a*(1-e^2)/(1+e*cos(theta0));   % km         % position vector length
h0=r0-Re;                         % km         % height 
V0=sqrt(mi*(2/r0-1/a));           % km/s       % orbital speed (ellipse)
V0 = 7.788; 

%de-orbiting retrograde burn
dV=0;                             % km/s       % impulsive delta V obtained
V0=V0-dV;                         % km/s       % effective inital speed 

%integration method
solv_kep='ode4';                  % solver method for orbital phase. ode4 is runge kutta
step_kep='30';  %s                % keplerian intergartor fixed step size   % 30 
stop_h=150;     %km               % thres alt for switch from Kep. to Atm. phase 

solv_atm = 'ode4';                % solver method for re-entry phase. ode4 is runge kutta
step_atm='0.1'; %s                % atmospheric integrator fixed step size     % 0.1

%conversions to I.S.
Re=Re*1000;                     %m
V0=V0*1000;                     %m/s
h0=h0*1000;                     %m
stop_h=stop_h*1000;             %m


h_0 = h0; 



%% Control conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
deg = 6.93; % Desired initial flight path angle.
% Apollo 4 initial flight path angle = 6.93 deg from horiz. % Taken Directly from the Apollo heating data for AS-501.
gamma_0 = deg*pi/180; %Convert to Radians

x_target = -4000000*10^3; % %Desired Final Range in km, AP4 final range 3797 km 

beta = 10^-11; %Step size--> 10^-11 used for ch. 6 tests.

n = 3; %Number of normalized time intervals. Apollo 4 = 19
Nmax = 1000; %Max number of iterations.
Imin = 1; %Max value of I for convergence.

delta_u(:,1) = -0.2*ones(n,1); %Initial step

t_f = 589.5; %Sensible guess for initial final time,


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Constraints
% % %
% Max heating and loading allowed, based on the max heating
% ... and loading from Apollo 4. Max loading experienced on Apollo 4 ... was roughly L_max = 7.25 g and max heat flux was q_max = 425 ... BTU/ft^2-s or 482.65 W/cm^2 or 4.83*10^6 W/m^2.
% ... See Sep 16th tech memo.
q_limit = 795.96*10^6; %max heat flux allowed in w/m^2 
L_limit = 12; %max loading allowed in g's.
h_limit = h_0; %Max altitude allowed in m.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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


%%
sys = 'controller_dynamics_OC';
kepOut=sim(sys,'Solver',solv_kep,'FixedStep',step_kep); %This runs the Simulink model using above values. 
%   updating initial conditions, 2nd run using ouput from 1st

%space environment
ndays1=kepOut.time(end)/60/60/24;  %days since intial time of simulation
jdate=jdate+ndays1;                %new julian date at beginning of 2nd run
F107_avg=90.85;         %SFU: updating
F107_day=86.0;          %SFU: updating
Kp=1;

%re-entry path initial values (at t=0 s)
x0=kepOut.x(end);            %m        %travelled distance
gamma0=kepOut.gamma(end);    %rad      %flight path angle
theta0=kepOut.theta(end);    %rad      %true anomaly
h0=kepOut.h(end);            %m        %height
V0=kepOut.V(end);            %m/s      %orbital speed


%Re-entry simulation, atmospheric phase, more accuracy is needed.
%Smaller integration step required (e.g. 0.1 s)
%Simulation stops when altitude reaches 0 km (default)
%%%%%%%%%%%%%%%%%%
stop_h=0;   %m          %final desired altitude
%%%%%%%%%%%%%%%%%%

atmOut = sim(sys,'Solver',solv_atm,'FixedStep',step_atm,'StartTime','kepOut.time(end)');


%% Results analysis and plotting
%   union of 1st and 2nd run results

h=[kepOut.h;atmOut.h];
x=[kepOut.x;atmOut.x];
time=[kepOut.time;atmOut.time];
V=[kepOut.V;atmOut.V];
Vx=[kepOut.Vx;atmOut.Vx];
Vz=[kepOut.Vz;atmOut.Vz];
gamma=[kepOut.gamma;atmOut.gamma];
theta=[kepOut.theta;atmOut.theta];
rho=[kepOut.rho;atmOut.rho];
dotV=[kepOut.dotV;atmOut.dotV];

%strutural loading [g]
L=abs(dotV/g0);             
% heat flux 
q_dot=1.83e-4*V.^3.*sqrt(rho./rcurv);

c = [h x time V Vx Vz gamma theta rho dotV L q_dot];
%   data filtering
%section to be revised
%excludes unsensed data due to discrete stopping criteria
%only (end) value is typically wrong (i.e. h(end)<0)


%%
n_sa=length(time);
k=0;
for i=n_sa:-1:1
    if h(i)<0 %|| h(i)>1.05*(r_a*1000-Re)
        k=k+1;
        h(i)=[];
        x(i)=[];
        time(i)=[];
        V(i)=[];
        Vz(i)=[];
        Vx(i)=[];
        gamma(i)=[];
        theta(i)=[];
        rho(i)=[];
        dotV(i)=[];
        q_dot(i)=[];
    end 
end
%%


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