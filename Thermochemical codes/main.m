%Run this file
%change alt1,alt2 to plot between those altitudes
load('Database/US1976_matlabdata.mat')
alt1 = 80000;
alt2 = 200000;
CdNewtonian = 4/3;

GammaM = zeros(2011,1);
Cpmax = zeros(2011,1);
Z = 79500;
while true
    Z = Z + 500;
    if Z > 1000000
        break
    end
    index = find(Zm==Z);
    T = TK(index);
    P = PPa(index); %in Pascals
    P = P / 100000; %converts it into BAR.

    gamma = getgamma(T,P);
    GammaM(index) = gamma;
    Cpmax(index) = ((gamma+1)^2 / (4*gamma))^(gamma*(gamma-1)) * (4 / (gamma+1));



end

Zmcrop = zeros(2011,1);
GammaMcrop = zeros(2011,1);
Cpmaxcrop = zeros(2011,1);

for index = 1:length(Zm)
    if alt2> Zm(index) && Zm(index) > alt1
    Zmcrop(index) = Zm(index);
    GammaMcrop(index) = GammaM(index);
    Cpmaxcrop(index) = Cpmax(index);
    end

end

Zmcrop = nonzeros(Zmcrop);
GammaMcrop = nonzeros(GammaMcrop);
Cpmaxcrop = nonzeros(Cpmaxcrop);


CdM = Cpmaxcrop .* CdNewtonian;


subplot(2,2,1)
plot(Zmcrop,GammaMcrop)
ylabel("Specific heat ratio of Air, γ")
xlabel("Height Z (m)")
title("γ vs Altitude, US1976 model with Bahadori numerical estimation of γ")


load("Database\SurfaceAreas.mat") %load theta surface area values.

subplot(2,2,2)

plot(ThetaM,AM/100)
ylabel("Normalized Surface Area normal to flow (cm2)")
xlabel("Angle (°)")
title("Projected area / Reference area (100) variation for orbits above 80km")






subplot(2,2,3)
plot(Zmcrop/1000,CdM)
ylabel("Drag coefficient")
xlabel("Height Z (km)")
title("Cd vs Altitude, Design2 (0° fin angle, modified Newtonian Theory, using γ plot)")


