function [T,P,rho,grav_const,MolecularWeight] = ATM1976(Z)
%ENTER ALTITUDE AS A MULTIPLE OF 500 IN METERS.
%OUTPUTS TEMPERATURE,PRESSURE,DENSITY in SI UNITS.
% addpath('/Users/zhanghanwen/Library/CloudStorage/OneDrive-Nexus365/3YP AtCUBESAT/Simulation/Gayhub/CubesatTrajectory/Thermochemical codes/Database')


load("Database/US1976_matlabdata.mat")

index = find(Zm==Z);

rho = rhokgm3(index);
T = TK(index);
P = PPa(index);
grav_const = g(index);
MolecularWeight = M(index);
