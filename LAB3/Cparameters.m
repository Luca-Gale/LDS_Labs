function out = Cparameters(in)
%CPARAMETERS computes the paramters of the compensator given the estimated
%model
%

%% Get the input variable
theta=in; % \hat theta_t

%% Modify the code ONLY here
% Nc=[3.04 -3.12 0.13]; % numerator of C(z)
% Dc =[1 -0.987 -0.013]; % denominator of C(z)

T_s = 0.001; % samplin time
t_r = 0.1; % rising time
O = 0.3; % overshoot

xi = abs(log(O)) / (sqrt(pi^2 + (log(O))^2 ));
omega_n = 1.8 / t_r;

p = exp( (  (-xi*omega_n) + 1j*omega_n*sqrt(1-xi^2)   )*T_s );
D_1 = [1 0 0];
D_2 = [1 -p];
D_3 = [1 -conj(p)];

D_star = conv(conv(D_1, D_2), D_3);

B = [theta(2), theta(3)];
A = conv([1 -2 1], [1, theta(1)]);

[D, Nc] = diophantine(A, B, D_star);
Dc = conv(D, [1,-1]);

%% Construct the outuput
out=[Nc Dc];


end






