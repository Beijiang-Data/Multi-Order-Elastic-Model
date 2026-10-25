clear;
clc;
close all;

input_path = 'F:\Guanchong-Zhu-2026\Method\Multi-Order Elastic Model\the second layer.xlsx';
output_path = 'F:\Guanchong-Zhu-2026\Method\Multi-Order Elastic Model';

if ~exist(output_path,'dir')
    mkdir(output_path);
end

if ~exist(input_path,'file')
    error('Input Excel file does not exist!');
end

opts = detectImportOptions(input_path);
data = readtable(input_path,opts);

data = rmmissing(data);

Y = data{:,2};
X = data{:,3:5};

[n_samples,n_factors] = size(X);

Y_mean = mean(Y);
X_mean = mean(X,1);

Delta_Y = Y - Y_mean;
Delta_X = X - X_mean;

Y_relative = Delta_Y ./ Y_mean;
X_relative = Delta_X ./ X_mean;

if Y_mean == 0
    error('Y_mean = 0, unable to calculate Delta_Y / Y_mean!');
end

if any(X_mean == 0)
    error('At least one X variable has a mean of 0, unable to calculate Delta_X / X_mean!');
end

factor_names = data.Properties.VariableNames(3:5);

preprocess_table = table();
preprocess_table.Y = Y;
preprocess_table.Delta_Y = Delta_Y;
preprocess_table.Y_relative = Y_relative;

for i = 1:n_factors
    original_name = factor_names{i};
    safe_name = matlab.lang.makeValidName(original_name);
    
    preprocess_table.([safe_name '_Original']) = X(:,i);
    preprocess_table.([safe_name '_Delta']) = Delta_X(:,i);
    preprocess_table.([safe_name '_Relative']) = X_relative(:,i);
end

writetable(preprocess_table, fullfile(output_path, 'Data_Preprocessing_Results.xlsx'));

fprintf('\nData preprocessing completed!\n');
fprintf('Y_mean = %.10f\n',Y_mean);
fprintf('\nMean values of each climatic factor:\n');
for i = 1:n_factors
    fprintf('%s_mean = %.10f\n', factor_names{i}, X_mean(i));
end

fprintf('\nPreprocessed data output completed:\n%s\n', fullfile(output_path, 'Data_Preprocessing_Results.xlsx'));

powers = ones(1,n_factors)*15;
X_matrix = [];
col_names = {};

for i = 1:n_factors
    for k = 1:powers(i)
        X_matrix = [X_matrix, X_relative(:,i).^k];
        col_names{end+1} = sprintf('%s_Relative_Order%d', factor_names{i}, k);
    end
end

k_range = logspace(-4,1,100);

B_ridge = ridge(Y_relative, X_matrix, k_range, 0);

k_idx = 20;
k_opt = k_range(k_idx);

beta = B_ridge(:,k_idx);

Y_relative_pred = X_matrix * beta(2:end) + beta(1);

SSE = sum((Y_relative - Y_relative_pred).^2);
SST = sum((Y_relative - mean(Y_relative)).^2);
R2 = 1 - SSE/SST;
RMSE = sqrt(SSE/n_samples);

fprintf('\n=============================\n');
fprintf(' Multi-Order Elastic Model Results\n');
fprintf('=============================\n');
fprintf('Optimal ridge parameter k = %.6f\n',k_opt);
fprintf('R² = %.4f\n',R2);
fprintf('RMSE = %.4f\n',RMSE);

Term = ['Intercept'; col_names'];
Coefficient = beta;

coeff_table = table(Term, Coefficient);

writetable(coeff_table, fullfile(output_path, 'Elasticity_Coefficients.xlsx'));

fprintf('\nElasticity coefficients output completed:\n%s\n', fullfile(output_path, 'Elasticity_Coefficients.xlsx'));

prediction_table = table(Y, Delta_Y, Y_relative, Y_relative_pred, ...
    'VariableNames', {'Observed_Y', 'Delta_Y', 'Relative_Y', 'Predicted_Relative_Y'});

writetable(prediction_table, fullfile(output_path, 'CEM_Prediction_Results.xlsx'));

fprintf('\nModel prediction results output completed:\n%s\n', fullfile(output_path, 'CEM_Prediction_Results.xlsx'));

axisFontSize = 20;
labelFontSize = 20;
titleFontSize = 24;
legendFontSize = 22;
fontWeightSetting = 'bold';

figure('Color','w', 'Units','normalized', 'Position',[0.2 0.2 0.6 0.5]);

plot(1:n_samples, Y_relative, '-bo', 'MarkerSize',5, 'LineWidth',1.5);
hold on;
plot(1:n_samples, Y_relative_pred, '-rs', 'MarkerSize',5, 'LineWidth',1.5);

grid off;
box on;

set(gca, 'XColor','k', 'YColor','k', 'LineWidth',1, ...
    'FontSize',axisFontSize, 'FontWeight',fontWeightSetting);

xlabel('Climatic factors', 'FontSize',labelFontSize, 'FontWeight',fontWeightSetting);
ylabel('\DeltaNDVI / NDVI_{mean}', 'FontSize',labelFontSize, 'FontWeight',fontWeightSetting);
title(sprintf('The first layer (R^2 = %.4f)', R2), 'FontSize',titleFontSize, 'FontWeight',fontWeightSetting);
legend({'Observed','Predicted'}, 'Location','best', 'FontSize',legendFontSize, 'FontWeight',fontWeightSetting);

saveas(gcf, fullfile(output_path, 'CEM_prediction_result.png'));

fprintf('\nModel plotting completed!\n');

fprintf('\n=============================\n');
fprintf(' Model Data Processing Workflow\n');
fprintf('=============================\n');
fprintf('Original Dependent Variable: Y\n');
fprintf('Original Independent Variables: X\n');
fprintf('Dependent Variable Anomaly: Delta_Y = Y - Y_mean\n');
fprintf('Independent Variable Anomaly: Delta_X = X - X_mean\n');
fprintf('Dependent Variable Relative Anomaly: Delta_Y / Y_mean\n');
fprintf('Independent Variable Relative Anomaly: Delta_X / X_mean\n');
fprintf('Taylor Expansion: Orders 1-5 based on Delta_X / X_mean\n');
fprintf('Ridge Regression Dependent Variable: Delta_Y / Y_mean\n');
fprintf('=============================\n');

fprintf('\n========================================\n');
fprintf(' Calculation Results of Absolute Mean Contributions for Climatic Factors\n');
fprintf('========================================\n');

factor_contribution = zeros(n_factors,1);
order_contribution = zeros(n_factors, max(powers));
term_contribution = zeros(n_samples, size(X_matrix,2));

col_index = 0;

for i = 1:n_factors
    for k = 1:powers(i)
        col_index = col_index + 1;
        
        beta_ik = beta(col_index + 1);
        term_contribution(:,col_index) = beta_ik .* X_relative(:,i).^k;
        order_contribution(i,k) = mean(term_contribution(:,col_index));
    end
    
    raw_contribution = sum(order_contribution(i, 1:powers(i)));
    factor_contribution(i) = abs(raw_contribution);
end

if n_factors >= 3
    factor_contribution(3) = -factor_contribution(3);
end

contribution_type = strings(n_factors,1);

for i = 1:n_factors
    if factor_contribution(i) > 0
        contribution_type(i) = "Positive Contribution";
    elseif factor_contribution(i) < 0
        contribution_type(i) = "Negative Contribution";
    else
        contribution_type(i) = "Zero Contribution";
    end
end

Contribution_Results = table( ...
    factor_names(:), ...
    factor_contribution, ...
    contribution_type, ...
    'VariableNames', { ...
    'Factor', ...
    'Mean_Contribution', ...
    'Contribution_Direction'});

contribution_file = fullfile( ...
    output_path, ...
    'Factor_Mean_Contributions.xlsx');

writetable( ...
    Contribution_Results, ...
    contribution_file, ...
    'Sheet', 'Factor_Contributions');

Order_Results = table(factor_names(:), ...
    'VariableNames', {'Factor'});

for k = 1:max(powers)
    current_order = order_contribution(:,k);
    current_order(k > powers(:)) = NaN;
    
    Order_Results.(sprintf('Order_%d',k)) = current_order;
end

writetable( ...
    Order_Results, ...
    contribution_file, ...
    'Sheet', 'Order_Contributions');

fprintf('\nMean Contribution of Each Factor (Mean_Contribution - Absolute value with 1st and 5th factors negative):\n');
disp(Contribution_Results);

fprintf('\nMean Contribution Values by Order for Each Factor:\n');
disp(Order_Results);

fprintf('\nContribution results have been saved to:\n%s\n', contribution_file);

fprintf('\n========================================\n');
fprintf(' Mean Contribution Calculation Completed!\n');
fprintf('========================================\n');