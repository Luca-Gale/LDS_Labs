clear
close all
clc


doublecheck = zeros(4,1);

for i=1:200

% Model generation / Data generation

N = 8000; % length of the data
Ts = 1;  % sampling time

% Coefficientes of the actual model (e.g. ARX model):
% numerators and denominators of the model with feedback
Fden = [1 -0.96 0.97];  
Fnum = [0 2.99 -0.2];    
Gden = [1 -0.96 0.97];
Gnum = [1 0 0];

% Input generation: e.g. sum of sinusoids whose frequencies are in the
% interval [pi*alpha1 pi*alpha2]
alpha1 = 1/8;
alpha2 = 1/2;
k_sin = 5; % number of sinusoids
u = idinput(N,'sine',[alpha1 alpha2],[],k_sin)*4;



% Generate NORMALIZED WGN noise en0 (i.e. en0 has variance equal to one)
en0 = idinput(N,'rgs');

% noise variance sigma0^2 of e0=sigma_0*en0
noiseVar = 4.6^2;

% Creation of the actual model (Ts is the sampling time, this parameter
% is optional and the default value is Ts = 1)
m0 = idpoly([],Fnum,Gnum,Gden,Fden,noiseVar,Ts); 

% Creation of the iddata object (a container for all the identification 
% data), in this case we store only the input data u
u = iddata([],u,Ts); 

% Creation of the iddata object, in this case we store only the normalized 
% noise data en0
en0 = iddata([],en0,Ts); 

% Generation of the output data given model, input and noise
y = sim(m0, [u en0]);

% Creation of an iddata object for the output and input data 
data = iddata(y,u);



%% PEM method

% Estimation of an ARX model

% orders of the ARX model
% orders_arx(1) = nA
% orders_arx(2) = nB
% orders_arx(3) = nk
orders_arx = [2 2 1]; 

% ARX model estimation
m_arx = arx(data,orders_arx);

% Plot the coefficients of the estimated model
m_arx.a % A(z)
m_arx.b % B(z)
m_arx.NoiseVariance % sigma^2


% Estimation of an ARMAX model

% orders of the ARMAX model
% orders_armax(1) = nA
% orders_armax(2) = nB
% orders_armax(3) = nC
% orders_armax(4) = nk
orders_armax = [2 2 1 1];

% ARMAX model estimation
m_armax = armax(data,orders_armax);

% Plot the coefficient of the estimated model
m_armax.a % A(z)
m_armax.b % B(z)
m_armax.c % C(z)
m_armax.NoiseVariance % sigma^2


% Estimation of an OE model

% orders of the OE model
% orders_oe(1) = nB
% orders_oe(2) = nF
% orders_oe(3) = nk
orders_oe = [2 1 1]; 

% OE model estimation
m_oe = oe(data,orders_oe);

% Plot the coefficient of the estimated model
m_oe.b % B(z)
m_oe.f % F(z)
m_oe.NoiseVariance % sigma^2


% Estimation of a Box-Jenkins model

% coefficients of the BJ model
% orders_bj(1) = nB
% orders_bj(2) = nC
% orders_bj(3) = nD
% orders_bj(4) = nF
% orders_bj(5) = nk
orders_bj = [2 1 1 2 1];

% BJ model generation
m_bj = bj(data,orders_bj);

% Plot the coefficient of the estimated model
m_bj.b % B(z)
m_bj.c % C(z)
m_bj.d % D(z) 
m_bj.f % F(z)
m_bj.NoiseVariance % sigma^2


%% Bode plot of the estimated input-output transfer functions
% Mettere label
% figure;
% bodeplot(m0,m_arx)
% figure;
% bodeplot(m0,m_armax)
% figure;
% bodeplot(m0,m_oe)
% figure;
% bodeplot(m0,m_bj)

%% Confidence interval
i = N;
S = zeros(4,4);

while i>0 
    S = S+sens(y.OutputData,u.InputData,i)*sens(y.OutputData,u.InputData,i)';
    i=i-1;
end

S = S./N;
P = m_arx.NoiseVariance*inv(S);

bounds = 1.96*[sqrt(P(1,1)/N); sqrt(P(2,2)/N); sqrt(P(3,3)/N); sqrt(P(4,4)/N)]

a = m_arx.a; % A(z)
b = m_arx.b; % B(z)

a10 = -0.96;
a20 = 0.97;
b00 = 2.99;
b10 = -0.2;

errors = [abs(a(2)-a10); abs(a(3)-a20); abs(b(2)-b00); abs(b(3)-b10)]

check = [errors(1)<=bounds(1); errors(2)<=bounds(2); errors(3)<=bounds(3); errors(4)<=bounds(4)]



doublecheck = doublecheck + check;

end


function Psi = sens(y,u,t)
    if t==2
        Psi = [-y(t-1); 0; u(t-1); 0];
    elseif t==1
        Psi = [0; 0; 0; 0];
    else
        Psi = [-y(t-1); -y(t-2); u(t-1); u(t-2)];
    end
end