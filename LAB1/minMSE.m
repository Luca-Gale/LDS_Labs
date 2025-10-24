function results = minMSE(lambdas, betas, y2, y1, x1, x2)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

N = size(lambdas,2);
results = zeros(N^2,3);
sigma_square_hat = 4.2;

count = 1;

for i=1:N
    for j=1:N
        beta = betas(j);
        lambda = lambdas(i);
        y_hat_2 = g_map(y1, x2, x1, beta, lambda, sigma_square_hat);
        mse = mean((y_hat_2 - y2).^2);
        results(count, 1) = beta;
        results(count, 2) = lambda;
        results(count, 3) = mse;
        count = count+1;
    end
end

end


function  y_hat= g_map(y1, x, x1, beta, lambda, sigma)
    y_hat = lambda*Cauchy_kernel(x, x1, beta)*inv(lambda*Cauchy_kernel(x1, x1, beta)+ ...
            sigma*eye(size(x1, 1)))*y1;
end