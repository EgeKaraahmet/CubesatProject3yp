%[rho0,h0,SH]=getSH(H) this function gives coefficient for an exponential
% atmosphere model
%works below 100 km
function [rho0,h0,SH]=getSH(H)
    %input data from table 8-4 of Fundamentals of Astrodynamics, D.A. Vallado
    h0=[0,25,30,40,50,60,70,80,90];
    rho0=[1.225,3.899e-2,1.774e-2,3.972e-3,1.057e-3,3.206e-4,8.770e-5,1.905e-5,3.396e-6];
    SH=[7.249,6.349,6.682,7.554,8.382,7.714,6.549,5.799,5.382];
    n=fix(H/10);
    if H>100 || H<0
        error('H must be in the range 0km < H < 100 km; current value H=%.2fkm\n',H) 
    end
    if H<30 
        n=2;
        if H<25 
            n=1;
        end 
    end
    if H==100
        n=9;
    end
    h0=h0(n);
    rho0=rho0(n);
    SH=SH(n);
end