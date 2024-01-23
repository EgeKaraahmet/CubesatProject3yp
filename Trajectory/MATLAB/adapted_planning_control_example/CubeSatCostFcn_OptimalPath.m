function J = CubeSatCostFcn_OptimalPath(stage, x, u)
   

  
    h = x(1);
    X = x(2);
    V = x(3);
    Re = 6371.0088 * 10^3;           %km         % Earth mean radius

    % Current position vector
    R = [X; h + Re];

    % Desired position vector at an altitude of h* = 30 km
    % h_star = 27.4188 * 10^3; % 30 km in meters --> A = 225 cm^2
    % X_star = -2551560;    %% seems not correct --> A = 225 cm^2

    h_star = 30 * 10^3;    %% @A=125cm^2
    X_star = -4.5871*10^9; %% interpolation from x_ref, h_ref @ A=125cm^2
    R_ref = [X_star; h_star + Re];
    q_ref = 600 * 10^3; 
    rcurv = 0.1; 

    % Weighting matrix P and Q
    P = [5   0;
         0   5];

    Q = 20; 

    % heat flux
    SH = 8397.5; 
    rho = 1.225 * exp(- h/(SH));    % my Python model, using SH = 8.43 and rho0=1.221 as constant
    q_max=1.83e-4*V.^3.*sqrt(rho/rcurv);

   
    % Cost function
    J = (R - R_ref)' * P * (R - R_ref) + (q_max - q_ref) * Q *(q_max-q_ref);
end

