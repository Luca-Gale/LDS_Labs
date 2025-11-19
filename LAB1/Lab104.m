clear all;
close all;
clc;

load(".\manipulator1.mat");

%% Data plotting

figure;
subplot(3,1,1);
plot(t, u1);
xlabel('Time (s)');
ylabel('Angular Position u1');
grid on;

subplot(3,1,2);
plot(t, u1d);
xlabel('Time (s)');
ylabel('Angular Velocity u1d');
grid on;

subplot(3,1,3);
plot(t, y1);
xlabel('Time (s)');
ylabel('Torque y1');
grid on;
sgtitle('First data set');

figure;
subplot(3,1,1);
plot(t, u2);
xlabel('Time (s)');
ylabel('Angular Position u2');
grid on;

subplot(3,1,2);
plot(t, u2d);
xlabel('Time (s)');
ylabel('Angular Velocity u2d');
grid on;

subplot(3,1,3);
plot(t, y2);
xlabel('Time (s)');
ylabel('Torque y2');
grid on;
sgtitle('Second data set');

%% Computation of the acceleration

u1dd = derivative(u1d, Ts);
u2dd = derivative(u2d, Ts);

%% Inverse dynamic

x1 = [u1 u1d u1dd];
x2 = [u2 u2d u2dd];

sigma_square = 4.2;
lambda = 10;
beta = 100;

%% Validation Step

y_hat_2 = g_map(y1, x2, x1, 100, 10, sigma_square);

figure;
hold on;
plot(t, y_hat_2, 'DisplayName', ...
    "\beta = " + 100 + ", \lambda = " + 10);

plot(t, y2, 'k', 'LineWidth', 1.2, 'DisplayName', 'y_2');

legend show;
xlabel('Time');
ylabel('Torque');
title('$y_2$ and $\hat{y}_2$', 'Interpreter', 'latex');

%% Validation with other values of the parameters

PARAM = {
    dictionary(["lambda","beta"], [10,  100]);
    dictionary(["lambda","beta"], [0.1, 100]);
    dictionary(["lambda","beta"], [10,  1])
    };

figure;
hold on;

for i = 1:numel(PARAM)
    values = PARAM{i};
    y_hat_2 = g_map(y1, x2, x1, values("beta"), values("lambda"), sigma_square);
    plot(t, y_hat_2, 'DisplayName', ...
        "\beta = " + values("beta") + ", \lambda = " + values("lambda"));
end

plot(t, y2, 'k', 'LineWidth', 1.2, 'DisplayName', 'y_2');

legend show;
xlabel('Time');
ylabel('Torque');
title('$y_2$ and $\hat{y}_2$', 'Interpreter', 'latex');

%% MSE comparison

MSE = zeros(size(PARAM, 1), 1);
for i = 1:numel(PARAM)
    values = PARAM{i};
    y_hat_2 = g_map(y1, x2, x1, values("beta"), values("lambda"), sigma_square);
    error = y_hat_2-y2;
    MSE(i) = (transpose(error)*error) / N;
end
MSE

%% Role of λ and β

% Role of Lambda
LAMBDA_COLLECTION = {
    dictionary(["lambda","beta"], [0.1,  100]);
    dictionary(["lambda","beta"], [1,    100]);
    dictionary(["lambda","beta"], [10,   100]);
    dictionary(["lambda","beta"], [100,  100]);
    dictionary(["lambda","beta"], [1000, 100]);
    };
figure;
hold on;
for i = 1:numel(LAMBDA_COLLECTION)
    values = LAMBDA_COLLECTION{i};
    y_hat_2 = g_map(y1, x2, x1, values("beta"), values("lambda"), sigma_square);
    plot(t, y_hat_2, 'DisplayName', ...
        "\beta = " + values("beta") + ", \lambda = " + values("lambda"));
end
plot(t, y2, 'k', 'LineWidth', 1.2, 'DisplayName', 'y_2');
legend show;
xlabel('Time');
ylabel('Torque');
title('$y_2$ and $\hat{y}_2$ changing $\lambda$', 'Interpreter', 'latex');

% Role of Beta
BETA_COLLECTION = {
    dictionary(["lambda","beta"], [100,    1]);
    dictionary(["lambda","beta"], [100,   10]);
    dictionary(["lambda","beta"], [100,  100]);
    dictionary(["lambda","beta"], [100, 1000]);
    };
figure;
hold on;
for i = 1:numel(BETA_COLLECTION)
    values = BETA_COLLECTION{i};
    y_hat_2 = g_map(y1, x2, x1, values("beta"), values("lambda"), sigma_square);
    plot(t, y_hat_2, 'DisplayName', ...
        "\beta = " + values("beta") + ", \lambda = " + values("lambda"));
end
plot(t, y2, 'k', 'LineWidth', 1.2, 'DisplayName', 'y_2');
legend show;
xlabel('Time');
ylabel('Torque');
title('$y_2$ and $\hat{y}_2$ changing $\beta$', 'Interpreter', 'latex');

%% Feedforward Control

w = [zeros(1,10) 10*ones(1,(N-10))]';
wd = derivative(w, Ts);
wdd = derivative(wd, Ts);

xf = [w wd wdd];

yf_matrix = zeros(size(w, 1), numel(PARAM));
for i = 1:numel(PARAM)
    values = PARAM{i};
    yf_matrix(:,i) = g_map(y1, xf, x1, values("beta"), values("lambda"), sigma_square);
end



function  y_hat = g_map(y1, x, x1, beta, lambda, sigma_square)
K11 = lambda * Cauchy_kernel(x1, x1, beta);
K21 = lambda * Cauchy_kernel(x, x1, beta);
y_hat = K21 * inv(K11 + sigma_square * eye(size(x1, 1))) * y1;
end

