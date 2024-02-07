
[rho_interp, SH_interp] = interpolate_exponential_model(10)

function [rho_interp, SH_interp] = interpolate_exponential_model(h_input)
   
    load("Database/US1976_matlabdata.mat")

    % Define the height values
    h0 = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240];

    % Initialize arrays to store results
    rho = zeros(size(h0));
    SH = zeros(size(h0));

    % Loop through each height value and calculate results
    for i = 1:length(h0)
        [rho(i), SH(i)] = exponential_model_data(h0(i));
    end

    % Interpolate for the given height values
    rho_interp = interp1(h0, rho, h_input, 'linear');
    SH_interp = interp1(h0, SH, h_input, 'linear');
    rho_0 = rho_interp/(exp((h_input)))
end
