% Define the height values
h = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240];

% Initialize arrays to store results
rho = zeros(size(h));
SH = zeros(size(h));

% Loop through each height value and calculate results
for i = 1:length(h)
    [rho(i), SH(i)] = exponential_model_data(h(i));
end

% Create a table
resultTable = table(h', rho', SH', 'VariableNames', {'Height', 'Density', 'SpecificHeat'});

% Display the table
disp(resultTable);
