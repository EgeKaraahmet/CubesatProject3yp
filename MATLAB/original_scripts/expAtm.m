%This function calculates density through an exponential model
%rho=expAtm(H)
%H is the geopotential height in km
%rho is the density in kg/m^3
%rho0,h0,SH are updated each 10 km (source: Wertz, 1978)
function rho=expAtm(H)
[rho0,h0,SH]=getSH(H);
rho=rho0*exp(-(H-h0)/SH);   %kg/m^3
end