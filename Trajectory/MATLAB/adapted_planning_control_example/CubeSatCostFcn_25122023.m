function J = CubeSatCostFcn_25122023(stage, x, u)
    % Extract relevant variables
    h = x(1);
    h_testing = evalin('base', 'h_atm_reference_signal');
    no_of_ratio = length(h_testing) / 12; 
    J = 0; 


    % Define the cost function
    for i = 1:12
        h_testing_tbe = h_testing(floor(i * no_of_ratio)); 
        J = J + (h_testing_tbe - h).^(2); 
    end

  

 

end


