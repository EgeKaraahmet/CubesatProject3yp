% This is a test by Ege to see if Matlab application can push changes to
% github. 5

% I want to show you the differences between cloud and your local files.

function state_derivative = motion_eqns(t, state)
    % Extract state variables
    x = state(1);
    y = state(2);
    z = state(3);
    vx = state(4);
    vy = state(5);
    vz = state(6);

    % constant
    R_Earth = 6371e3;     % Radius of the Earth (km)
    GM = 3.986004418e14;  % m^3 / s^2
    earth = Earth(R_Earth, GM); 
    database = MSISE_90(); % solar activity counted data

    % Earth's rotation angular velocity (rad / s)
    omega = [0;0;7.2921159e-5];  

    % Speed of sound 
    a = 343;    % m/s
   

    %% Compute vectors v and r based on time t
    v = [vx; vy; vz];
    r = [x; y; z];

    %% Compute density at a given altitude 
    rho = earth.density(r);

    %% Absolute velocity 
    V = v - cross(omega,[x;y;z]);
    V_mag = norm(V); 
    
    %% Ballistic coefficient
    M = 3   ;   % mass of the Cubesat
    S = 20  ;   % Total surface area! need to consider aerobrake

    mach_no = V_mag/a; 
    [~,Cd]=database.get_mach_data(mach_no); 

    Cb = Cd * S / M;  % ! This is not correct 

 

    %% Derivatives 
    % Compute derivatives
    dvdt = -GM * r / norm(r)^3 - 0.5 * rho * Cb * norm(v)^2 * v;
    drdt = V;

    % Create the state derivative vector
    state_derivative = [drdt ; dvdt];
end


