load("A127.mat")
x_127 = abs(x);
h_127 = h; 
time_127 = time/60; 

load("A126p5.mat")
x_126p5 = abs(x);
h_126p5 = h; 
time_126p5 = time/60; 

load("A126.mat")
x_126 = abs(x);
h_126 = h; 
time_126 = time/60;

load("A125p5.mat")
x_125p5 = abs(x);
h_125p5 = h; 
time_125p5 = time/60;

load("A125.mat")
x_125 = abs(x);
h_125 = h; 
time_125 = time/60; 

load("A124p5.mat")
x_124p5 = abs(x);
h_124p5 = h; 
time_124p5 = time/60; 

load("A124.mat")
x_124 = abs(x);
h_124 = h; 
time_124 = time/60;

load("A123p5.mat")
x_123p5 = abs(x);
h_123p5 = h; 
time_123p5 = time/60;

load("A123.mat")
x_123 = abs(x);
h_123 = h; 
time_123 = time/60; 

figure
plot(x_127,h_127, ...
    x_126p5,h_126p5, ...
    x_126,h_126, ...
    x_125p5,h_125p5, ...
    x_125,h_125, ...
    x_124p5,h_124p5, ...
    x_124,h_124, ...
    x_123p5,h_123p5, ...
    x_123,h_123)
ylabel('Altitude (km)')
xlabel('Horizontal distance travelled (km)')

figure
plot(time_127,h_127, ...
    time_126p5,h_126p5, ...
    time_126,h_126, ...
    time_125p5,h_125p5, ...
    time_125,h_125, ...
    time_124p5,h_124p5, ...
    time_124,h_124, ...
    time_123p5,h_123p5, ...
    time_123,h_123)
ylabel('Altitude (km)')
xlabel('Time (hours)')

%% Calculations 
altitude = 6378; 

x_125_landing = x_125(end);
x_125_landing = x_125_landing - floor(x_125_landing/(2*pi*altitude));
position_ref= sqrt(altitude^2+x_125_landing^2);


% small perturbation
x_127_landing = x_127(end);
x_127_landing = x_127_landing - floor(x_127_landing/(2*pi*altitude))
position_127= sqrt(altitude^2+x_127_landing^2);
pe_127 = abs(position_127-position_ref);

x_126_landing = x_126(end);
x_126_landing = x_126_landing - floor(x_126_landing/(2*pi*altitude));
position_126= sqrt(altitude^2+x_126_landing^2);
pe_126 = abs(position_126-position_ref);

x_124_landing = x_124(end);
x_124_landing = x_124_landing - floor(x_124_landing/(2*pi*altitude));
position_124= sqrt(altitude^2+x_124_landing^2);
pe_124 = abs(position_124-position_ref);

x_123_landing = x_123(end);
x_123_landing = x_123_landing - floor(x_123_landing/(2*pi*altitude));
position_123= sqrt(altitude^2+x_123_landing^2);
pe_123 = abs(position_123-position_ref);

x_126p5_landing = x_126p5(end);
x_126p5_landing = x_126p5_landing - floor(x_126p5_landing/(2*pi*altitude));
position_126p5= sqrt(altitude^2+x_126p5_landing^2);
pe_126p5 = abs(position_126p5-position_ref);

x_125p5_landing = x_125p5(end);
x_125p5_landing = x_125p5_landing - floor(x_125p5_landing/(2*pi*altitude));
position_125p5= sqrt(altitude^2+x_125p5_landing^2);
pe_125p5 = abs(position_125p5-position_ref);

x_124p5_landing = x_124p5(end);
x124p5_landing = x_124p5_landing - floor(x_124p5_landing/(2*pi*altitude));
position_124p5= sqrt(altitude^2+x_124p5_landing^2);
pe_124p5 = abs(position_124p5-position_ref);

x_123p5_landing = x_123p5(end);
x_123p5_landing = x_123p5_landing - floor(x_123p5_landing/(2*pi*altitude));
position_123p5= sqrt(altitude^2+x_123p5_landing^2);
pe_123p5 = abs(position_123p5-position_ref);

% Define positional error values
positional_errors = [pe_123, pe_123p5, pe_124, pe_124p5, pe_125p5, pe_126, pe_126p5, pe_127];

% Create a table
T_pe = table(positional_errors', 'VariableNames', {'Positional_Error_km'}, 'RowNames', {'123', '123.5', '124', '124.5', '125.5', '126', '126.5', '127'});

% Display the table
disp(T_pe)
