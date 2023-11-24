%This function combines Jacchia J71 and Exponential Atmosphere for %calculating density in the band [0;2500] km
%A transition function should be introduced to avoid discontinuities
function rho = unirho(h,jdate,F107_avg,F107_act,Kp) 
     h=h/1000; %conversion [m] to [km]
    rho=0;
    %from 0 to 100 km of altitude use Exponential Atmosphere
    if h<100 && h>=0
        rho=expAtm(h);
    end
    %from 100 to 2500 km of altitude use Jacchia J71
    if h>=100 && h<=2500 
       rho=J71_density_simulink(h,jdate,F107_avg,F107_act,Kp);
    end
    %A warning is printed if a negative altitude is predicted by %the simulation (due to Simulink discrete stopping criterion)
    if h<0
       fprintf('Warning: mismatched density, altitude: %3.2e km',h) 
       rho=1.225;
    end
end