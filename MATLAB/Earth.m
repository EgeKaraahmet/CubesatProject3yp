classdef Earth
    % Base class for the Earth, which allows us to 
    %       1) do coordinate system conversions 
    %       2) distance dependent gravitation acceleration is assumed
    %       3) the USA Standard Atmosphere model was assumed: 
 

    properties
        radius
        grav_k
    end

    methods
        function obj = Earth(radius, grav_k)
            % R_earth = 6371 km
            obj.radius = radius;

            % grav_k = GM 
            obj.grav_k = grav_k;
        end

        %% ECI: polar <-> cartesian 
        function [x,y,z] = cartesian(obj, lat, lon, alt)
            % Convert spherical coordinates (latitude, longitude, altitude) to 3D Cartesian coordinates (x, y, z).
            r = alt + obj.radius;
            x = r * cos(lat) * cos(lon);
            y = r * cos(lat) * sin(lon);
            z = r * sin(lat);
        end

        function [lat,lon,alt] = polar(obj, x, y, z)
            % Convert 3D Cartesian coordinates (x, y, z) to spherical coordinates (latitude, longitude, altitude).
            r = norm([x, y, z]);
            lat = asin(z / r);
            lon = atan2(y, x);
            alt = r - obj.radius;
        end


        %% gravity acceleartion in the ECI frame 
        function g = gravity(obj, r)  
            g = -(obj.grav_k / (r^2.0));
        end

        %% the USA Standard Atmosphere model
        function alt = altitude(obj, r)
            alt = norm(r) - obj.radius;
        end

        function rho = density(obj,r)
            alt = norm(r) - obj.radius;
            rho = 1.221 * exp(-alt / 8.43e3);
        end

        

        
    end     
end

