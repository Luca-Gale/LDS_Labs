clear all
close all
clc

%% Model generation / Data generation

N = 1000; % length of the data
Ts = 1;  %  sampling time

% load data
load house.mat;

%% Data Visualization
% Define time axis based on the number of samples N
t = 1:N; 

figure(1);
plot(t, u, 'b', t, y, 'r');
title('Input (u) and Output (y) Comparison');
xlabel('Time [samples]');
ylabel('Amplitude');
legend('Input (u)', 'Output (y)');
grid on;

figure(2);
plot(t, u, 'b');
title('Input Signal (u)');
xlabel('Time [samples]');
ylabel('Amplitude');
legend('u');
grid on;

figure(3);
plot(t, y, 'r');
title('Output Signal (y)');
xlabel('Time [samples]');
ylabel('Amplitude');
legend('y');
grid on;

%% detrending
% Creation of an iddata object for the output and input data
data = iddata(y,u);

% Compute the means of the data and store them in mu.InputOffset and
% mu.OutputOffset
mu = getTrend(data,0);

% Perform the data detrend
data_d = detrend(data,mu);

%% Create models
na = 50;
nb = 50;
nk = 1;

% M1
Orders = [na nb nk];
Option = arxRegulOptions('RegulKernel', 'DI');
[Lambda, R, log_negative_likelihood_m1, df_m1] = arxRegul(data_d, Orders, Option); % use data_d
arxOpt_m1 = arxOptions;
arxOpt_m1.Regularization.Lambda = Lambda;
arxOpt_m1.Regularization.R = R;
model_m1 = arx(data_d, Orders, arxOpt_m1);

% M2
Orders = [na nb nk];
Option = arxRegulOptions('RegulKernel', 'TC');
[Lambda, R, log_negative_likelihood_m2, df_m2] = arxRegul(data_d, Orders, Option); % use data_d
arxOpt_m2 = arxOptions;
arxOpt_m2.Regularization.Lambda = Lambda;
arxOpt_m2.Regularization.R = R;
model_m2 = arx(data_d, Orders, arxOpt_m2);

% M3
Orders = [na nb nk];
Option = arxRegulOptions('RegulKernel', 'SS');
[Lambda, R, log_negative_likelihood_m3, df_m3] = arxRegul(data_d, Orders, Option); % use data_d
arxOpt_m3 = arxOptions;
arxOpt_m3.Regularization.Lambda = Lambda;
arxOpt_m3.Regularization.R = R;
model_m3 = arx(data_d, Orders, arxOpt_m3);

% M4
Orders = [na nb nk];
m = 5;
[Lambda, R, log_negative_likelihood_m4, df_m4] = arxRegulRF2(data_d , Orders, m); % use data_d
arxOpt_m4 = arxOptions;
arxOpt_m4.Regularization.Lambda = Lambda;
arxOpt_m4.Regularization.R = R;
model_m4 = arx(data_d, Orders, arxOpt_m4);

%% Impulse Response Visualization (A(z) and B(z))
models = {model_m1, model_m2, model_m3, model_m4};
model_names = {'M1 (DI)', 'M2 (TC)', 'M3 (SS)', 'M4 (RF)'};

for i = 1:4
    figure('Name', sprintf('Impulse Responses for %s', model_names{i}), 'NumberTitle', 'off');
    
    % Plot Impulse Response of A(z)
    subplot(2, 1, 1);
    plot(models{i}.A, 'LineWidth', 1.5);
    title([model_names{i}, ': Impulse Response of A(z)']);
    xlabel('Lag');
    ylabel('Amplitude');
    grid on;
    
    % Plot Impulse Response of B(z)
    subplot(2, 1, 2);
    plot(models{i}.B, 'LineWidth', 1.5, 'Color', 'r');
    title([model_names{i}, ': Impulse Response of B(z)']);
    xlabel('Lag');
    ylabel('Amplitude');
    grid on;
end

%% Compare different negative log likelihood
negative_log_likelihood = [log_negative_likelihood_m1;
                            log_negative_likelihood_m2;
                            log_negative_likelihood_m3;
                            log_negative_likelihood_m4;]

%% Hold cross validation
input_training = u(1:500);
output_training = y(1:500);

input_validation = u(501:1000);
output_validation = y(501:1000);

% Creation of an iddata object for the output and input data 
data_training = iddata(output_training, input_training);
data_validation = iddata(output_validation, input_validation);

% Compute the means of the data and store them in mu.InputOffset and
% mu.OutputOffset
mu_training = getTrend(data_training, 0);
mu_validation = getTrend(data_validation, 0);

% Perform the data detrend
data_training_d = detrend(data_training, mu_training);
data_validation_d = detrend(data_validation, mu_validation);

% ARMAX model estimation
m1_training = arx(data_training_d, Orders, arxOpt_m1);
m2_training = arx(data_training_d, Orders, arxOpt_m2);
m3_training = arx(data_training_d, Orders, arxOpt_m3);
m4_training = arx(data_training_d, Orders, arxOpt_m4);

opt = compareOptions('InitialCondition','z');

models_training = {m1_training, m2_training, m3_training, m4_training};
model_labels    = {'M1 (DI)', 'M2 (TC)', 'M3 (SS)', 'M4 (RF)'};
horizons        = 1:4;

for i = 1:length(models_training)
    figure('Name', sprintf('Validation: %s', model_labels{i}), 'NumberTitle', 'off');
    
    for h = horizons
        subplot(2, 2, h);
        
        compare(data_validation_d, models_training{i}, h, opt);
        
        title(sprintf('%s: %d-step ahead', model_labels{i}, h));
        grid on;
    end
end

%% Cross correlation and auto correlation
% the confidence interval corresponds to 99%
models = {model_m1, model_m2, model_m3, model_m4};
model_names = {'M1 (DI)', 'M2 (TC)', 'M3 (SS)', 'M4 (RF)'};
figure('Name', 'Residual Analysis Comparison', 'NumberTitle', 'off');

for i = 1:4
    subplot(2, 2, i);
    resid(models{i}, data_d, 'corr'); 
    
    title(['Residuals: ', model_names{i}]);
    grid on;
end

%% SURE index
SURE_indexes = zeros(4, 1);
SURE_indexes(1) = surek(model_m1, df_m1);
SURE_indexes(2) = surek(model_m2, df_m2);
SURE_indexes(3) = surek(model_m3, df_m3);
SURE_indexes(4) = surek(model_m4, df_m4);
SURE_indexes

%% computation time
comp_times = zeros(20,2);
Orders = [na nb nk];
Option = arxRegulOptions('RegulKernel', 'SS');
m = 5;

for i=1:20
    tic

    % M3
    [Lambda, R, log_negative_likelihood_m3, df_m3] = arxRegul(data_d, Orders, Option); % use data_d

    comp_times(i,1) = toc;
end

std_SS = std(comp_times(:,1));
mean_SS = mean(comp_times(:,1));

for i=1:20
    tic

    % M4
    [Lambda, R, log_negative_likelihood_m4, df_m4] = arxRegulRF2(data_d , Orders, m); % use data_d

    comp_times(i,2) = toc;
end

std_RF = std(comp_times(:,2));
mean_RF = mean(comp_times(:,2));

%% Results table
Method = {'M3 (SS)'; 'M4 (RF)'};
MeanTime = [mean_SS; mean_RF];
StdTime  = [std_SS; std_RF];

ResultsTable = table(MeanTime, StdTime, ...
    'RowNames', Method, ...
    'VariableNames', {'Mean', 'Std'});

disp(ResultsTable)
