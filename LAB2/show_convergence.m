clear
close all
clc

Ts = 1;

% Actual model
Fden = [1 -0.96 0.97];
Fnum = [0 2.99 -0.2];
Gden = [1 -0.96 0.97];
Gnum = [1 0 0];

alpha1 = 1/8;
alpha2 = 1/2;
k_sin = 5; % number of sinusoids

% NOise variance
noiseVar = 4.6^2;

% Create actual model
m0 = idpoly([],Fnum,Gnum,Gden,Fden,noiseVar,Ts);

% orders of the ARX model
orders_arx = [2 2 1];

% orders of the ARMAX model
orders_armax = [2 2 1 1];

% orders of the OE model
orders_oe = [2 1 1];

% coefficients of the BJ model
orders_bj = [2 1 1 2 1];

structures(4) = struct('name', [], 'orders', [], 'function', []);

% Initialize the structures array with model names and orders
structures(1) = struct('name', 'ARX', 'orders', orders_arx, 'function', @arx);
structures(2) = struct('name', 'ARMAX', 'orders', orders_armax, 'function', @armax);
structures(3) = struct('name', 'OE', 'orders', orders_oe, 'function', @oe);
structures(4) = struct('name', 'BJ', 'orders', orders_bj, 'function', @bj);

for j=1:4

    % Initialise figure for bodeplot
    figure;
    h = bodeplot(m0);
    opts = getoptions(h);
    hold on;

    % Make m0 thicker
    kids = findall(h, 'type', 'line');
    set(kids, 'LineWidth', 2.5);

    for N = round(linspace(200, 8000, 10))

        % Input for actual model
        u = idinput(N,'sine',[alpha1 alpha2],[],k_sin)*4;
        u = iddata([],u,Ts);

        % Normalised wgn
        en0 = idinput(N,'rgs');
        en0 = iddata([],en0,Ts);

        % Simulate output
        y = sim(m0, [u en0]);

        % Creation of an iddata object for the output and input data
        data = iddata(y,u);

        % ARX model estimation
        m = structures(j).function(data, structures(j).orders);
        bodeplot(m, opts);
    end

end