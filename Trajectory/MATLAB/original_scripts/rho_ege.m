function [rho, SH] = rho_ege(H)
    % Define the data points
    heights = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, ...
               100, 110, 120, 130, 140, 150, 160, 170, ...
               180, 190, 200, 210];  % Corrected range
    rhos = [1.225, 1.2238, 1.2227, 1.2215, 1.2204, 1.2192, ...
            1.2181, 1.2169, 1.2158, 1.2146, 1.2135, 1.2123, ...
            1.2112, 1.21, 1.2088, 1.2077, 1.2065, 1.2054, ...
            1.2042, 1.2031, 1.2019, 1.2008]; % Adjusted the number of elements
    SHs = [8434.1, 8432.5, 8430.6, 8428.7, 8426.8, 8425, ...
           8423.1, 8421.2, 8419.3, 8417.5, 8415.6, 8413.7, ...
           8411.8, 8410, 8408.1, 8406.2, 8404.3, 8402.4, ...
           8400.6, 8398.7, 8396.8, 8394.9]; % Adjusted the number of elements

    % Perform interpolation
    rho = interp1(heights, rhos, H, 'spline');
    SH = interp1(heights, SHs, H, 'spline');
    
    % Check if H is within the valid range
    if H < 0 || H > 210
        error('H must be in the range 0 km < H < 210 km; current value H=%.2f km\n', H)
    end
end




% function rho=rho_ege(H)
%     %input data from 
%     if  H<0
%         error('H must be in the range 0km < H < 100 km; current value H=%.2fkm\n',H) 
%     end
%     if H>=0  && H<10
%         rho0 = 1.225;
%         SH = 8434.1;
%     end
% 
%     if H>=10  && H<20
%         rho = 1.2238;
%         SH = 8432.5;
%     end
% 
%     if H>=20  && H<30
%         rho = 1.2227;
%         SH = 8430.6;
%     end
% 
%     if H>=30  && H<40
%         rho = 1.2215;
%         SH = 8428.7;
%     end
% 
%     if H>=40  && H<50
%         rho = 1.2204;
%         SH = 8426.8;
%     end
% 
%     if H>=50  && H<60
%         rho = 1.2192;
%         SH = 8425;
%     end
% 
%     if H>=60  && H<70
%         rho = 1.2181;
%         SH = 8423.1;
%     end
% 
%     if H>=70  && H<80
%         rho = 1.2169;
%         SH = 8421.2;
%     end
% 
%     if H>=80  && H<90
%         rho = 1.2158;
%         SH = 8419.3;
% 
%     end
% 
%     if H>=90  && H<100
%         rho = 1.2146;
%         SH = 8417.5;
%     end
%     if H>=100  && H<110
%         rho = 1.2135;
%         SH = 8415.6;
%     end
%     if H>=110  && H<120
%         rho = 1.2123;
%         SH = 8413.7;
%     end
% 
%     if H>=120  && H<130
%        rho = 1.2112;
%         SH = 8411.8;
%     end
% 
%     if H>=130  && H<140
%         rho = 1.21;
%         SH = 8410;
%     end
% 
%     if H>=140  && H<150
%        rho = 1.2088;
%         SH = 8408.1;
%     end
% 
%     if H>=150  && H<160
%        rho = 1.2077;
%         SH = 8406.2;
%     end
% 
%     if H>=160  && H<170
%        rho = 1.2065;
%         SH = 8404.3;
%     end
% 
%     if H>=170  && H<180
%        rho = 1.2054;
%         SH = 8402.4;
%     end
% 
%     if H>=180  && H<170
%        rho = 1.2042;
%         SH = 8400.6;
%     end
% 
%     if H>=190  && H<200
%         rho = 1.2031;
%         SH = 8398.7;
%     end
% 
%     if H>=200  && H<210
%         rho = 1.2019;
%         SH = 8396.8;
%     end
% 
% 
% 
% end