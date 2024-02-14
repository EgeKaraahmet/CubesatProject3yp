
load("A125.mat")
x_125 = x;
h_125 = h; 
time_125 = time; 

load("A124p5.mat")
x_124p5 = x;
h_124p5 = h; 
time_124p5 = time; 

load("A124.mat")
x_124 = x;
h_124 = h; 
time_124 = time;

load("A123p5.mat")
x_123p5 = x;
h_123p5 = h; 
time_123p5 = time;

load("A123.mat")
x_123 = x;
h_123 = h; 
time_123 = time; 

plot(x_125,h_125, ...
    x_124p5,h_124p5, ...
    x_124,h_124, ...
    x_123p5,h_123p5, ...
    x_123,h_123)