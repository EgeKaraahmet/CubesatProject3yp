%This function computes the density using Jacchia J71 model
%rho=J71_density_simulink(H,jdate,F107_avg,F107_act,Kp)
%input: height [90-2500 km], julian date, F10.7 3-monthly average, F10.7
%three hourly average, Kp planetary index
%output: density value
%Reference: Satellite Orbits, section 3.5.3
function rho=J71_density_simulink(H,jdate,F107_avg,F107_act,Kp)

    %% Exosperic Temperature 
    T_c = 379.0+3.24*F107_avg+1.3*(F107_act-F107_avg);
    % Geographic corrections
    % The actual exospheric temperature depends on the local hour angle
    % of the Sun with respect to the satellite.
    % It also depends, however, on the declination of the Sun and the
    % geographic latitude of the satellite.
    % The actual exospheric temperature T1 with the diurnal variations
    % included can be computed from a more compliated formula available
    % in the referenced book.
    % Geomagnetic corrections
    % using the three-hourly planetary geomagnetic index Kp for a time
    % 6.7 hours earlier than the time under consideration
    f = 0.5*(tanh(0.04*(H-350)+1));
    dT_infH = 28.0*Kp+0.03*exp(Kp);
    dT_infL = 14.0*Kp+0.02*exp(Kp);
    dT_gm = f*dT_infH+(1-f)*dT_infL;
    %transition function
    %H>350km
    %H<350km
    T_e = T_c+dT_gm;   %geomagnetic corrected exospheric temperature
    
    %% Standard Density 
    %Use bi-polynomial fit from Gill (1996)
    coeff=getCoeff(H,T_e); %coefficient matrix
    logRho=0;
    for i=0:5
        for j=0:4
            logRho = logRho + coeff(i+1,j+1)*((H/1000)^i)*((T_e/1000)^j);
        end 
    end

    %% Corrections 
    %Correction due to semi-annual density variation in thermosphere
    Phi=(jdate - 2400000.5 -36204)/365.2422; %number of tropical years since Jan1, 1958
    
    tauSA=Phi + 0.09544*((0.5+0.5*sin(2*pi*Phi+6.035))^1.65 - 0.5);
    fZ=(5.876e-7*H^2.331 + 0.06328)*exp(H*-2.868e-3);
    gt=0.02835 + 0.3817*(1+0.4671*sin(2*pi*tauSA+4.137))*sin(4*pi*tauSA+4.259);
    logRhoSA=fZ*gt;
    %Correction due to geomagnetic activities
     if H<350
            logRhoGM = (0.012*Kp + 1.2e-5*exp(Kp))*(1-f);%Transition function from temp calc
     else
        logRhoGM=0;
     end
    %Seasonal-latitude correction
    %logRhoSL=0.014*(Z-90).*exp(-0.0013*(Z-90).^2)*sin(1*pi*Phi -1.72)*sin(lat)^3/abs(sin(lat));
    logRho=logRho+logRhoGM+logRhoSA;
    rho=10^logRho;
end
