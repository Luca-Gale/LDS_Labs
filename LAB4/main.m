clear all;
close all;
clc;

%% Model generation / Data generation

N = 1000; % length of the data
Ts = 1;  %  sampling time

% load data
load house.mat;

%% Data plotting
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

%% Data detrending

% Creation of an iddata object for the output and input data
data = iddata(y, u, Ts);

% Compute the means of the data and store them in mu.InputOffset and
% mu.OutputOffset
mu = getTrend(data, 0);

% Perform the data detrend
data_d = detrend(data, mu);

% Delay estimation
nk = 1;

%% Model structure design

% M1 ARMAX
% orders of the ARMAX model
orders_armax_1(1) = 1;
orders_armax_1(2) = 1;
orders_armax_1(3) = 1;
orders_armax_1(4) = nk;

% ARMAX model estimation
m_armax_1 = armax(data_d, orders_armax_1);

% M2 ARMAX
% orders of the ARMAX model
orders_armax_2(1) = 2;
orders_armax_2(2) = 2;
orders_armax_2(3) = 2;
orders_armax_2(4) = nk;

% ARMAX model estimation
m_armax_2 = armax(data_d, orders_armax_2);

% M3 OE
% orders of the OE model
orders_oe(1) = 2;
orders_oe(2) = 2;
orders_oe(3) = nk;

% OE model estimation
m_oe = oe(data_d, orders_oe);

% M4 BJ
% coefficients of the BJ model
orders_bj(1) = 2;
orders_bj(2) = 2;
orders_bj(3) = 2;
orders_bj(4) = 2;
orders_bj(5) = nk;

% BJ model generation
m_bj = bj(data_d, orders_bj);

%% Residual Analysis

% the confidence interval corresponds to 99%
figure;
resid(m_armax_1, data_d, 'corr', 'r');
figure;
resid(m_armax_2, data_d, 'corr', 'r');
figure;
resid(m_oe, data_d, 'corr', 'r');
figure;
resid(m_bj, data_d, 'corr', 'r');

%% Zero-Pole cancellation analysis

% Zeros and Poles plot
figure;
iopzplot(m_armax_1)
figure;
iopzplot(m_armax_2)
figure;
iopzplot(m_oe)
figure;
iopzplot(m_bj)

%% Criteria with complexity term

% Model names
Model = {'ARMAX 1'; 'ARMAX 2'; 'OE'; 'BJ'};

% SURE index
SURE = [
    sure(m_armax_1)
    sure(m_armax_2)
    sure(m_oe)
    sure(m_bj)
    ];

% AIC index
AIC = [
    aic(m_armax_1)
    aic(m_armax_2)
    aic(m_oe)
    aic(m_bj)
    ];

% BIC index
BIC = [
    bic(m_armax_1)
    bic(m_armax_2)
    bic(m_oe)
    bic(m_bj)
    ];

ResultsTable = table(Model, SURE, AIC, BIC)


%% Hold-out Cross-validation

% prediction k-steps ahead from zero initial conditions
k = 1;

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
m_armax_training_1 = armax(data_training_d, orders_armax_1);
m_armax_training_2 = armax(data_training_d, orders_armax_2);
m_oe_training = oe(data_training_d, orders_oe);
m_bj_training = bj(data_training_d, orders_bj);

opt = compareOptions('InitialCondition','z');

figure;
compare(data_validation_d, m_armax_training_1, k, opt)
figure;
compare(data_validation_d, m_armax_training_2, k, opt)
figure;
compare(data_validation_d, m_oe_training, k, opt)
figure;
compare(data_validation_d, m_bj_training, k, opt)

