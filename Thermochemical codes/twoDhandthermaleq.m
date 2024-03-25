clc
clear all

T5(1) = input("Enter T5: ");
T4(1) = input("Enter T4: ");
T3(1) = input("Enter T3: ");
T2(1) = input("Enter T2: ");

error = 10;

k = 237;
sigma = 5.670373*10^-8;
e = 0.2;
rho = 2700;
c = 890;

A5 = 2.5*10^-4;
A4 = 1.6575*10^-3; %whole area
A3 = 1.9*10^-3; %whole area
A2 = 0.01;


m5 = rho * A5;
m4 = rho * A4;
m3 = rho * A3;
m2 = rho * A2;


L5 = 0.105;
L4 = 0.195;
L3 = 0.201;
L2 = 0.198;

t = 0; %time is 0, in seconds
%1 orbit is assumed to be 90 minutes, which is roughly a true number.
dt = 0.1;
i = 1;
j = 0; %this counter counts how many times t reaches 2700.
orbittype = "Shadow";

while j < 80
   if  mod(j,2) == 0
       orbittype = "Sun";
   else
       orbittype = "Shadow";
   end
   if orbittype == "Shadow"
        q5in = 22.15;
        q3in = 22.15; %one part of q3
        q2in = 3.98;
   elseif orbittype == "Sun"
       q5in = 27;
       q3in = 33.892; %one part of q3
       q2in = 3.98;
   else
       stop("Error: Orbittype not defined")
   end

q5out = e * sigma * (T5(i)^4) * L5;
q4out = e * sigma * (T4(i)^4) * L4;
q3out = e * sigma * (T3(i)^4) * L3; %one part of q3
q2out = e * sigma * (T2(i)^4) * L2;

q32 = (k * (T3(i)-T2(i)) * 0.0095) / 0.1; %heat transfer from one fin of region 3
q24 = (k * (T2(i)-T4(i)) * 0.0085 * 2) / 0.05; %region 2 to region 4 whole.
q45 = (k * (T4(i)-T5(i)) * 0.0085 * 2) / 0.05; %region 4 to region 5 whole.




q5net = q5in - q5out + q45;
q4net = -q4out - q45 + q24;
q3net = (q3in - q3out - q32) * 2;   %this is all two parts of region 3.
q2net = (q32*2) - q2out - q24 + q2in;


T5(i+1) = T5(i) + (q5net*dt/(m5*c));
T4(i+1) = T4(i) + (q4net*dt/(m4*c));
T3(i+1) = T3(i) + (q3net*dt/(m3*c));
T2(i+1) = T2(i) + (q2net*dt/(m2*c));
time(i+1) = t + j*2700;
t = t+dt;

if t >= 2700
    j = j + 1;
    t = t - 2700;
end

i = i+1;

end

timeminutes = time / 60;

plot(timeminutes,T2,timeminutes,T3,timeminutes,T4,timeminutes,T5)
legend('Region 2','Region 3','Region 4','Region 5')
xlabel('Time (min)')
ylabel('Temperature (K)')
title('Temperature variation during multiple orbits from an initial position of 290K')

