classdef MSISE_90
%%% This database output CL and Cd, 
% The MSIS model is to describe the Earth's atmosphere and 
% its behavior, particularly in the upper atmosphere and ionosphere

% The "MSISE-90" model is used 

    properties (Access = private)
        mach_data
        atm_data
    end
    
    methods
        function obj = MSISE_90()
            % Mach number data
            obj.mach_data = [
                0 12.6 0.0;
                0.2 12.6 0.2;
                0.5 12.6 0.3;
                0.6 13.0 0.4;
                0.8 13.8 0.43;
                1.05 17.5 0.47;
                1.5 14.7 0.55;
                2.0 13.8 0.51;
                2.5 11.4 0.49;
                3.0 10.1 0.47;
                4.0 9.3 0.43;
                5.0 7.5 0.39;
                6.0 6.8 0.39;
                7.0 6.8 0.39;
                8.0 6.8 0.39;
                9.0 6.8 0.39;
                10.0 6.8 0.39;
            ];
            % Atm data 
            obj.atm_data = [
                0 300.2511 1.17E+00 1.01E+05 28.9502;
                20 206.2085 9.49E-02 5.62E+03 28.9502;
                40 257.6979 4.07E-03 3.02E+02 28.9502;
                60 244.1212 3.31E-04 2.32E+01 28.9502;
                80 196.3636 1.68E-05 9.45E-01 29.0175;
                100 184.0160 5.08E-07 2.81E-02 27.7137;
                120 374.9715 1.80E-08 2.17E-03 25.8745;
                140 635.5703 3.26E-09 7.03E-04 24.5349;
                160 787.5532 1.18E-09 3.31E-04 23.4225;
                180 877.6729 5.51E-10 1.80E-04 22.4106;
                200 931.2806 2.91E-10 1.05E-04 21.4734;
                220 963.2701 1.66E-10 6.44E-05 20.6108;
                240 982.4191 9.91E-11 4.09E-05 19.8292;
                260 993.9173 6.16E-11 2.66E-05 19.1337;
                280 1000.8427 3.94E-11 1.77E-05 18.5256;
                300 1005.0267 2.58E-11 1.20E-05 18.0015;
                320 1007.5620 1.72E-11 8.20E-06 17.5537;
                340 1009.1030 1.16E-11 5.69E-06 17.1721;
                360 1010.0423 7.99E-12 3.98E-06 16.8449;
                380 1010.6166 5.55E-12 2.81E-06 16.5597;
                400 1010.9688 3.89E-12 2.01E-06 16.3044;
                420 1011.1853 2.75E-12 1.44E-06 16.0669;
                440 1011.3190 1.96E-12 1.04E-06 15.8360;
                460 1011.4014 1.40E-12 7.55E-07 15.6008;
                480 1011.4526 1.01E-12 5.53E-07 15.3508;
                500 1011.4845 7.30E-13 4.07E-07 15.0760;
                520 1011.5043 5.31E-13 3.03E-07 14.7669;
                540 1011.5168 3.88E-13 2.27E-07 14.4148;
                560 1011.5245 2.85E-13 1.71E-07 14.0125;
                580 1011.5294 2.11E-13 1.31E-07 13.5547;
                600 1011.5325 1.56E-13 1.01E-07 13.0389;
                620 1011.5345 1.17E-13 7.89E-08 12.4665;
                640 1011.5357 8.79E-14 6.24E-08 11.8428;
                660 1011.5365 6.65E-14 5.01E-08 11.1779;
                680 1011.5370 5.08E-14 4.07E-08 10.4854;
                700 1011.5374 3.91E-14 3.36E-08 9.7818;
                720 1011.5375 3.04E-14 2.82E-08 9.0847;
                740 1011.5377 2.39E-14 2.39E-08 8.4111;
                760 1011.5377 1.90E-14 2.06E-08 7.7753;
                780 1011.5378 1.53E-14 1.79E-08 7.1884;
                800 1011.5378 1.25E-14 1.58E-08 6.6572;
                820 1011.5378 1.03E-14 1.40E-08 6.1849;
                840 1011.5379 8.64E-15 1.26E-08 5.7711;
                860 1011.5379 7.32E-15 1.14E-08 5.4132;
                880 1011.5379 6.28E-15 1.04E-08 5.1066;
                900 1011.5379 5.46E-15 9.47E-09 4.8460;
            ];
        end

        function [t, d, p, m] = get_atmospheric_data(obj, alt)
            altitude = obj.atm_data(:, 1);
            temperature = obj.atm_data(:, 2);
            density = obj.atm_data(:, 3);
            pressure = obj.atm_data(:, 4);
            molwt = obj.atm_data(:, 5);
            
            t = interpolate(altitude, temperature, alt);
            d = interpolate(altitude, density, alt);
            p = interpolate(altitude, pressure, alt);
            m = interpolate(altitude, molwt, alt);
        end

        function [cl, cd] = get_mach_data(obj, mach_no)
            ar_mach = obj.mach_data(:, 1);
            ar_cl = obj.mach_data(:, 2);
            ar_cd = obj.mach_data(:, 3);

            cl = interpolate(ar_mach, ar_cl, mach_no);
            cd = interpolate(ar_mach, ar_cd, mach_no);
        end
    end
end

function result = interpolate(x_list, y_list, x)
    if any(x_list(2:end) <= x_list(1:end-1))
        error('x_list must be in strictly ascending order!');
    end
    
    i = find(x_list <= x, 1, 'last');
    
    if x <= x_list(1)
        result = y_list(1);
    elseif x >= x_list(end)
        result = y_list(end);
    else
        result = y_list(i) + (x - x_list(i)) * (y_list(i + 1) - y_list(i)) / (x_list(i + 1) - x_list(i));
    end
end



