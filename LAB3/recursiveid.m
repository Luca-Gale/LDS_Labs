    function out=recursiveid(in)
%RECURSIVEID computes the recursive PEM estimate using the ARX model 
%            structure with na=1,nb=2, nk=0
%
%

%% Get the input variables
na=1;
nb=2;
theta=in(1:na+nb); % \hat theta_t-1 
ypast=in(na+nb+1:na+2*nb-1); % y(t-1)
upast=in(na+2*nb:2*na+2*nb-1); % u(t-1)
S=reshape(in(2*(na+nb):2*(na+nb)+(na+nb)^2-1),na+nb,na+nb); % S_t-1
u=in(2*(na+nb)+(na+nb)^2); % u(t)
y=in(2*(na+nb)+1+(na+nb)^2); % y(t)

%% Modify the code ONLY here
%Sn = eye(3);
% thetan=[-0.2 0.02 0.001]'; %\hat theta_t

lambda = .999;
phi = [-ypast, u, upast]';
Sn = lambda*S+phi*phi';
thetan = theta + inv(Sn)*phi*(y - phi'*theta);

%% Construct the outuput
out=[thetan; y; u; reshape(Sn,(na+nb)^2,1)]; 
