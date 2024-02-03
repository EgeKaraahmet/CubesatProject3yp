function [rho,SH] = exp_model(Z)
%ENTER ALTITUDE AS A MULTIPLE OF 500 IN METERS.
%OUTPUTS TEMPERATURE,PRESSURE,DENSITY in SI UNITS.
addpath('/Users/zhanghanwen/Library/CloudStorage/OneDrive-Nexus365/3YP CUBESAT/Simulation/Gayhub/CubesatTrajectory/Thermochemical codes/Database')
load("US1976_matlabdata.mat")



if mod(Z,500) ~= 0
    remainder = mod(Z,500);
    Z1 = Z - remainder;
    Z2 = Z1 + 500;
    index1 = find(Zm==Z1);
    index2 = find(Zm==Z2);
    rho2 = rhokgm3(index2);
    T2 = TK(index2);
    P2 = PPa(index2);
    grav_const2 = g(index2);
    MolecularWeight2 = M(index2);
else
    Z1 = Z;
    remainder = 0;
    index1 = find(Zm==Z1);
end



weight = remainder / 500;





rho1 = rhokgm3(index1);
T1 = TK(index1);
P1 = PPa(index1);
grav_const1 = g(index1);
MolecularWeight1 = M(index1);




if mod(Z,500) ~= 0
    rho = ((1-weight)*rho1)+((weight)*rho2);
    T = ((1-weight)*T1+(weight)*T2);
    P = ((1-weight)*P1+(weight)*P2);
    grav_const = ((1-weight)*grav_const1 + (weight)*grav_const2);
    MolecularWeight = ((1-weight)*MolecularWeight1 + (weight)*MolecularWeight2);
else
    rho = rhokgm3(index1);
    T = TK(index1);
    P = PPa(index1);
    grav_const = g(index1);
    MolecularWeight = M(index1);
end



R = 8.314 / (10^(-3)); 
SH = R.*T./(MolecularWeight.*grav_const);
