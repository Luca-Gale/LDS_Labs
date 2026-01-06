clear all
close all
clc

%% Model generation / Data generation

N = 1000; % length of the data
Ts = 1;  %  sampling time

% load data
load house.mat;

% %% plotting
% t = 1:N;
% figure(1);
% hold on;
% plot(t, u);
% plot(t, y);
% legend('u', 'y');
% 
% figure(2);
% plot(t, u);
% legend('u');
% 
% figure(3);
% plot(t, y);
% legend('y');

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
nk = 1;


%% M1 ARMAX
% Estimation of an ARMAX model

% orders of the ARMAX model
orders_armax_1(1) = 1;
orders_armax_1(2) = 1;
orders_armax_1(3) = 1;
orders_armax_1(4) = nk;


% ARMAX model estimation
m_armax_1 = armax(data_d,orders_armax_1);

% Plot the coefficient of the estimated model
m_armax_1.a; % A(z)
m_armax_1.b; % B(z)
m_armax_1.c; % C(z)
m_armax_1.NoiseVariance; % sigma^2

%% M2 ARMAX
% Estimation of an ARMAX model

% orders of the ARMAX model
orders_armax_2(1) = 2;
orders_armax_2(2) = 2;
orders_armax_2(3) = 2;
orders_armax_2(4) = nk;


% ARMAX model estimation
m_armax_2 = armax(data_d,orders_armax_2);

% Plot the coefficient of the estimated model
m_armax_2.a; % A(z)
m_armax_2.b; % B(z)
m_armax_2.c; % C(z)
m_armax_2.NoiseVariance; % sigma^2


%% Output error
% Estimation of an OE model

% orders of the OE model
orders_oe(1) = 2;
orders_oe(2) = 2;
orders_oe(3) = nk;


% OE model estimation
m_oe = oe(data_d,orders_oe);

% Plot the coefficient of the estimated model
m_oe.b; % B(z)
m_oe.f; % F(z)
m_oe.NoiseVariance; % sigma^2

%% Box Jenkins
% Estimation of a Box-Jenkins model

% coefficients of the BJ model
orders_bj(1) = 2;
orders_bj(2) = 2;
orders_bj(3) = 2;
orders_bj(4) = 2;
orders_bj(5) = nk;


% BJ model generation
m_bj = bj(data_d,orders_bj);

% Plot the coefficient of the estimated model
m_bj.b; % B(z)
m_bj.c; % C(z)
m_bj.d;% D(z) 
m_bj.f; % F(z)
m_bj.NoiseVariance; % sigma^2


%% Model structure determination

% Residual Analysis

% the confidence interval corresponds to 99%
% figure;
% resid(m_armax_1, data_d,'corr', 'r');
% figure;
% resid(m_armax_2, data_d,'corr', 'r');
% figure;
% resid(m_oe, data_d,'corr', 'r');
% figure;
% resid(m_bj, data_d,'corr', 'r');




% Zero-Pole cancellation analysis

% % Zeros and Poles plot  
% figure;
% iopzplot(m_armax_1)
% figure;
% iopzplot(m_armax_2)
% figure;
% iopzplot(m_oe)
% figure;
% iopzplot(m_bj)


% SURE index
SURE_indexes = zeros(4, 1);
SURE_indexes(1) = sure(m_armax_1);
SURE_indexes(2) = sure(m_armax_2);
SURE_indexes(3) = sure(m_oe);
SURE_indexes(4) = sure(m_bj);
SURE_indexes;

% AIC index
AIC_indexes = zeros(4, 1);
AIC_indexes(1) = aic(m_armax_1);
AIC_indexes(2) = aic(m_armax_1);
AIC_indexes(3) = aic(m_oe);
AIC_indexes(4) = aic(m_bj);
AIC_indexes;

% BIC index
BIC_indexes = zeros(4, 1);
BIC_indexes(1) = bic(m_armax_1);
BIC_indexes(2) = bic(m_armax_2);
BIC_indexes(3) = bic(m_oe);
BIC_indexes(4) = bic(m_bj);
BIC_indexes;

% Hold out cross validation
% (ATTENTION: for the cross-validation the comparison has to be done on 
% the validation dataset)
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
mu_training = getTrend(data_training,0);
mu_validation = getTrend(data_validation, 0);

% Perform the data detrend
data_training_d = detrend(data_training, mu_training);
data_validation_d = detrend(data_validation, mu_validation);
% temp_1 = data_d.InputData(1:500);
% temp_2 = data_d.InputData(501:1000);
% temp_3 = data_d.OutputData(1:500);
% temp_4 = data_d.OutputData(501:1000);
% data_training_d = iddata(temp_3, temp_1);
% data_validation_d = iddata(temp_4, temp_2);

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

