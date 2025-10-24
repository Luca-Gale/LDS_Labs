clear all;
close all;
clc;

load(".\manipulator1.mat");

%% Plot data

figure;
subplot(3,1,1)
plot(t, u1);
xlabel('Time (s)');
ylabel('Angular Position u1');
grid on;

subplot(3,1,2)
plot(t, u1d);
xlabel('Time (s)');
ylabel('Angular Velocity u1d');
grid on;

subplot(3,1,3)
plot(t, y1);
xlabel('Time (s)');
ylabel('Torque y1');
grid on;

figure;
subplot(3,1,1)
plot(t, u2);
xlabel('Time (s)');
ylabel('Angular Position u2');
grid on;

subplot(3,1,2)
plot(t, u2d);
xlabel('Time (s)');
ylabel('Angular Velocity u2d');
grid on;

subplot(3,1,3)
plot(t, y2);
xlabel('Time (s)');
ylabel('Torque y2');
grid on;

%% Computation of the acceleration

ud1_dot = derivative(u1d, Ts);
ud2_dot = derivative(u2d, Ts);

%% Inverse dynamic

x1 = [u1 u1d ud1_dot];
x2 = [u2 u2d ud2_dot];

sigma_square_hat = 4.2;
lambda = 10;
beta = 100;

%% Validation Step
param = {
    dictionary(["lambda","beta"], [10,  100]);
    dictionary(["lambda","beta"], [0.1, 100]);
    dictionary(["lambda","beta"], [10,  1])
};

figure; hold on;

for i = 1:numel(param)
    values = param{i}; 
    y_hat_2 = g_map(y1, x2, x1, values("beta"), values("lambda"), sigma_square_hat);
    plot(t, y_hat_2, 'DisplayName', ...
        "β = " + values("beta") + ", λ = " + values("lambda"));
end

plot(t, y2, 'k', 'LineWidth', 1.2, 'DisplayName', 'Real torque');
legend show;
xlabel('Time');
ylabel('Torque');
title('Validation Step');

%% Optimal value of hyperparameters (not required)
lambdas = round(logspace(0, 6, 50));
betas = round(logspace(0, 6, 50));
res = minMSE(lambdas, betas, y2, y1, x1, x2);

[~, idx] = min(res(:,3));
beta_opt = res(idx, 1);
lambda_opt = res(idx, 2); 
minMSE = res(idx, 3);
% Display optimal hyperparameters
fprintf('Optimal lambda: %.2f\n', lambda_opt);
fprintf('Optimal beta: %.2f\n', beta_opt);
fprintf('Minimum MSE: %.2f\n', minMSE);

figure
y_hat_2 = g_map(y1, x2, x1, beta_opt, lambda_opt, sigma_square_hat);
plot(t, y_hat_2, 'DisplayName', ...
        "β = " + beta_opt + ", λ = " + lambda_opt);
hold on;
plot(t, y2, 'k', 'LineWidth', 1.2, 'DisplayName', 'Real torque');
legend show;
xlabel('Time');
ylabel('Torque');
title('Optimal Parameters');

%%
wt = transpose([zeros(1,10) 10*ones(1,(N-10))]);
wtd = derivative(wt, Ts);
wtdd = derivative(wtd, Ts);

xf = [wt wtd wtdd];

lambda = 10;
beta = 100;

yf = g_map(y1, xf, x1, beta, lambda, sigma_square_hat);


function  y_hat= g_map(y1, x, x1, beta, lambda, sigma)
    y_hat = lambda*Cauchy_kernel(x, x1, beta)*inv(lambda*Cauchy_kernel(x1, x1, beta)+ ...
            sigma*eye(size(x1, 1)))*y1;
end
