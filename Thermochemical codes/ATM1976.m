function [T,P,rho] = ATM1976(Z)
%ENTER ALTITUDE AS A MULTIPLE OF 500 IN METERS.
%OUTPUTS TEMPERATURE,PRESSURE,DENSITY in SI UNITS.

load("US1976_matlabdata.mat")

index = find(Zm==Z);

rho = rhokgm3(index);
T = TK(index);
P = PPa(index);
