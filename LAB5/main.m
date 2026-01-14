clear all
close all
clc

%% Model generation / Data generation

N = 1000; % length of the data
Ts = 1;  %  sampling time

% load data
load house.mat;

%% plotting
t = 1:N;
figure(1);
hold on;
plot(t, u);
plot(t, y);
legend('u', 'y');

figure(2);
plot(t, u);
legend('u');

figure(3);
plot(t, y);
legend('y');

%% detrending

% Creation of an iddata object for the output and input data
data = iddata(y,u);

% Compute the means of the data and store them in mu.InputOffset and
% mu.OutputOffset
mu = getTrend(data,0);

% Perform the data detrend
data_d = detrend(data,mu);

% Delay estimation
%nk = delayest(data_d);

%% Create models
na = 50;
nb = 50;
nk = 1;

% M1
Orders = [na nb nk]; % [n_a, n_b, n_k]
Option = arxRegulOptions('RegulKernel', 'DI');
[Lambda, R, log_negative_likelihood_m1, df_m1] = arxRegul(data_d, Orders, Option); % use data_d
arxOpt_m1 = arxOptions;
arxOpt_m1.Regularization.Lambda = Lambda;
arxOpt_m1.Regularization.R = R;
model_m1 = arx(data_d, Orders, arxOpt_m1);

% M2
Orders = [na nb nk]; % [n_a, n_b, n_k]
Option = arxRegulOptions('RegulKernel', 'TC');
[Lambda, R, log_negative_likelihood_m2, df_m2] = arxRegul(data_d, Orders, Option); % use data_d
arxOpt_m2 = arxOptions;
arxOpt_m2.Regularization.Lambda = Lambda;
arxOpt_m2.Regularization.R = R;
model_m2 = arx(data_d, Orders, arxOpt_m2);

% M3
Orders = [na nb nk]; % [n_a, n_b, n_k]
Option = arxRegulOptions('RegulKernel', 'SS');
[Lambda, R, log_negative_likelihood_m3, df_m3] = arxRegul(data_d, Orders, Option); % use data_d
arxOpt_m3 = arxOptions;
arxOpt_m3.Regularization.Lambda = Lambda;
arxOpt_m3.Regularization.R = R;
model_m3 = arx(data_d, Orders, arxOpt_m3);

% M4
Orders = [na nb nk]; % [n_a, n_b, n_k]
m = 5;
[Lambda, R, log_negative_likelihood_m4, df_m4] = arxRegulRF2(data_d , Orders, m); % use data_d
arxOpt_m4 = arxOptions;
arxOpt_m4.Regularization.Lambda = Lambda;
arxOpt_m4.Regularization.R = R;
model_m4 = arx(data_d, Orders, arxOpt_m4);

%% Plot A(z), B(z)

% m1
figure;
plot(model_m1.A)
figure;
plot(model_m1.B)

% m2
figure;
plot(model_m2.A)
figure;
plot(model_m2.B)

% m3
figure;
plot(model_m3.A)
figure;
plot(model_m3.B)

% m4
figure;
plot(model_m4.A)
figure;
plot(model_m4.B)

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
mu_training = getTrend(data_training,0);
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

% figure;
% compare(data_validation_d, m1_training, 1, opt)
% figure;
% compare(data_validation_d, m1_training, 2, opt)
% figure;
% compare(data_validation_d, m1_training, 3, opt)
% figure;
% compare(data_validation_d, m1_training, 4, opt)
% figure;
% compare(data_validation_d, m2_training, 1, opt)
% figure;
% compare(data_validation_d, m2_training, 2, opt)
% figure;
% compare(data_validation_d, m2_training, 3, opt)
% figure;
% compare(data_validation_d, m2_training, 4, opt)
% figure;
% compare(data_validation_d, m3_training, 1, opt)
% figure;
% compare(data_validation_d, m3_training, 2, opt)
% figure;
% compare(data_validation_d, m3_training, 3, opt)
% figure;
% compare(data_validation_d, m3_training, 4, opt)
% figure;
% compare(data_validation_d, m4_training, 1, opt)
% figure;
% compare(data_validation_d, m4_training, 2, opt)
% figure;
% compare(data_validation_d, m4_training, 3, opt)
% figure;
% compare(data_validation_d, m4_training, 4, opt)
% 

%% Cross correlation and auto correlation
% the confidence interval corresponds to 99%
% figure;
% resid(model_m1, data_d,'corr', 'r');
% figure;
% resid(model_m2, data_d,'corr', 'r');
% figure;
% resid(model_m3, data_d,'corr', 'r');
% figure;
% resid(model_m4, data_d,'corr', 'r');

%% SURE index
% SURE index
SURE_indexes = zeros(4, 1);
SURE_indexes(1) = surek(model_m1, df_m1);
SURE_indexes(2) = surek(model_m2, df_m2);
SURE_indexes(3) = surek(model_m3, df_m3);
SURE_indexes(4) = surek(model_m4, df_m4);
SURE_indexes;

%% computation time
comp_time = zeros(21,1);
Orders = [na nb nk]; % [n_a, n_b, n_k]
Option = arxRegulOptions('RegulKernel', 'SS');
m = 5;

for i=1:21
    tic

    % M3
    [Lambda, R, log_negative_likelihood_m3, df_m3] = arxRegul(data_d, Orders, Option); % use data_d

    % M4
    [Lambda, R, log_negative_likelihood_m4, df_m4] = arxRegulRF2(data_d , Orders, m); % use data_d

    comp_time(i) = toc;
end

standard_deviation = std(comp_time)
mean_time = mean(comp_time)

