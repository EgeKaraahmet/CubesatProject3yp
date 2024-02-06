%this is a linear average value for temperatures,,altitude,g,molecular
%weight

load('Database/US1976_matlabdata.mat');

index1 = find(Zm==200000);
index2 = find(Zm==80000);




AvgTK = (TK(index1) + TK(index2)) / 2;
Avgg = (g(index1) + g(index2)) / 2;
AvgM = (M(index1) + M(index2)) / 2;

R = 8.314 / (1*10^-3);

SH = (R * AvgTK) / (AvgM * Avgg)