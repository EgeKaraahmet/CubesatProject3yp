function y_value = findCorrespondingY(x, y, target_x)
    % Check if the input vectors have the same length
    if length(x) ~= length(y)
        error('Input vectors x and y must have the same length');
    end
    
    % Find the index of the target_x in vector x
    index = find(x == target_x, 1);
    
    % If target_x is not found, perform linear interpolation
    if isempty(index)
        lower_index = find(x < target_x, 1, 'last');
        upper_index = find(x > target_x, 1, 'first');
        
        % Check if interpolation is possible
        if isempty(lower_index) || isempty(upper_index)
            error('Target value not found in vector x, and interpolation is not possible');
        end
        
        % Linear interpolation
        x_lower = x(lower_index);
        x_upper = x(upper_index);
        y_lower = y(lower_index);
        y_upper = y(upper_index);
        
        % Interpolate to find the corresponding value in vector y
        slope = (y_upper - y_lower) / (x_upper - x_lower);
        y_value = y_lower + slope * (target_x - x_lower);
    else
        % Use the index to find the corresponding value in vector y
        y_value = y(index);
    end
end

