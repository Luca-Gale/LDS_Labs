function  [Lambda, R, obj, df] = arxRegulRF2(data,Orders,mred)
%ARXREGULRF Determine the SS regularization matrix for an ARX model 
%           using random features.
%
%  [Lambda, R, obj] = arxRegulb(data,Orders,Opt)
%   
%  Variables Lambda, R, obj, data, Orders are the same of arxRegul.m
%
%  Mattia Zorzi December 2019

warning off
%% input data
u = data.u;
y = data.y;
if sum(Orders<0)>0 
    error('Orders should be a vector with nonnegative values')
end
if sum(Orders(1:2)==0)==2
    error('n_a and n_b cannot be zero together')
end
if ne(Orders(1),Orders(2)) & Orders(1)>0
    error('This case has not been implemented')
end
if Orders(1)==Orders(2) 
    lab='yy';
    p=Orders(1);
end
if Orders(1)==0
    lab='ny';
    p=Orders(2);
end
if Orders(2)==0
    lab='yn';
    p=Orders(1);
end
if ne(Orders(3),1)
    error('Arxregulb is only implemented with n_k=1')
end
 

%% kernel tipe (TC)
D=-eye(p);
for k=1:p-1;
    D(k+1,k)=1;
end
D=D*D;
Di=(D^-1)';


%% construct regression matrix
[Phi,m,nv,y]=Abuild(y,u,p,lab);

%% new output
switch lab
    case 'ny'
        Phin=Phi*Di;
        yn = Phin'*y;
    case 'yn'
        Phin=Phi*Di;
        yn = Phin'*y;
    case 'yy'
        Phin(:,:,1)=Phi(:,:,1)*Di;
        Phin(:,:,2)=Phi(:,:,2)*Di;
        yn(:,1) = Phin(:,:,1)'*y;
        yn(:,2) = Phin(:,:,2)'*y;
end
%% random features


%% optimization & random features
switch lab
    case 'ny'
        lambda_init=[1 1];
        lb=[0 0];
        Z=randn(p,mred);
    case 'yn'
        lambda_init=[1 1];
        lb=[0 0];
        Z=randn(p,mred);
    case 'yy'
        lambda_init=[1 1 1 1];
        lb=[0 0 0 0];
        Z(:,:,1)=randn(p,mred);
        Z(:,:,2)=randn(p,mred);
end
opt=optimoptions(@fmincon,'Display','off');
[lambda, obj] = fmincon(@(lambda) loglikb(lambda,Phin,yn,nv,Z,p,lab),lambda_init,[],[],[],[],lb,[],[],opt);
                                          
%% kernel 

switch lab
    case 'ny'
        B=Di*kernel(lambda(1),lambda(2),p)*Z;
        K=B*B'+10^-10*eye(size(B,1));
    case 'yn'
        B=Di*kernel(lambda(1),lambda(2),p)*Z;
        K=B*B'+10^-10*eye(size(B,1));
    case 'yy'
        B1=Di*kernel(lambda(1),lambda(3),p)*Z(:,:,1);
        B2=Di*kernel(lambda(2),lambda(4),p)*Z(:,:,2);
        K=zeros(2*p,2*p);
        K(1:p,1:p)=B1*B1'+10^-10*eye(size(B1,1));
        K(p+1:2*p,p+1:2*p)=B2*B2'+10^-10*eye(size(B2,1));
end

K=K/nv;
LK=chol(K)';
LKi=inv(LK);
R=LKi'*LKi;
Lambda=norm(R);
R=R/norm(R);

switch lab
    case 'ny'
        obj = obj+0.5*(nv^-1*y'*y+(size(y,1)-mred)*log(nv));
        df=trace(Phi*(Phi'*Phi+nv*K^-1)^-1*Phi');
    case 'yn'
        obj = obj+0.5*(nv^-1*y'*y+(size(y,1)-mred)*log(nv));
        df=trace(Phi*(Phi'*Phi+nv*K^-1)^-1*Phi');
    case 'yy'
        obj = obj+0.5*(nv^-1*y'*y+(size(y,1)-2*mred)*log(nv));
        df=trace([Phi(:,:,1) Phi(:,:,2)]*([Phi(:,:,1)'*Phi(:,:,1) Phi(:,:,1)'*Phi(:,:,2); 
            Phi(:,:,2)'*Phi(:,:,1) Phi(:,:,2)'*Phi(:,:,2)]+nv*K^-1)^-1*[Phi(:,:,1)'; Phi(:,:,2)']);
end


 
function [A,m,nv,ystacked]=Abuild(y,u,p,lab)
% ABUILD Compute regressor and noise variance of the model
%
% A(z)y(t)=B(z)u(t-1)+e(t)
%
% where A(z)=1+sum_{k=1}^na a_kz^-k and B(z)=sum_{k=0}^nb-1 b_kz^-k.
% NB the output is reduced to initialize the one step ahead predictor 
%
% INPUT
% y   output data 
% u   input data 
% p   =na=nb
% lab choose the model structure: 'yy' BJ NP model
%                                 'ny' OE NP model 
%                                 'yn' time series
%
%OUTPUT
%
% A     regression matrix with 3 indexes:
%          A(:,:,1)-A(:,:,my) regressors associated to outputs 1-my
%          A(:,:,my+1)-A(:,:,my+mu) regressors associated to inputs 1-mu
% m     third dimension of A
% nv    noise variance (different in each output channel)
% yred  reduced output with outputs stacked




%n = (dimension of y) - p, data in y^+

fat=2/3;%establishes ARX order
ARXorder=p;%maximum ARX order



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% computing A,T,m,yred, ytrue %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n=size(y,1)-p;
my=size(y,2);%number of outputs
mu=size(u,2);%number of inputs
T=[];
if lab=='yy'%ARX
    m=my+mu;
    for i=1:my;
        for j=1:n;
            f=y(:,i)';
            A(j,:,i)=fliplr(f([j:(j+p-1)]));
        end
        T=[T A(:,:,i)];
    end
    for i=(my+1):m;
        for j=1:n;
            f=u(:,i-my)';
            A(j,:,i)=fliplr(f([j:(j+p-1)]));
        end
        T=[T A(:,:,i)];
    end
elseif lab=='ny';%output error
    m=mu;
    for i=1:mu;
        for j=1:n;
            f=u(:,i)';
            A(j,:,i)=fliplr(f([j:(j+p-1)]));
        end
        T=[T A(:,:,i)];
    end
elseif lab=='yn';%time series
    m=my;
    for i=1:my;
        for j=1:n;
            f=y(:,i)';
            A(j,:,i)=fliplr(f([j:(j+p-1)]));
        end
        T=[T A(:,:,i)];
    end
end
yred = y(p+1:n+p,:);

ystacked = reshape(yred,n*my,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% computing nv (innovation variance) %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pr=min(ARXorder,floor((n*fat)/m));
pr=min(p,pr);
Tr=[];
for i=1:m;
    Tr=[Tr A(:,1:pr,i)];
end
for i=1:my;
    [X] = lscov(Tr,yred(:,i));
    SD=sqrt(sum((yred(:,i)-Tr*X).^2)/(n-size(Tr,2)));
    nv(i)=SD^2;
end

%q=pem([y u]);nv=q.NoiseVariance;

% in the case that the noise variance is the same for each output channel
%nv=mean(nv);
%SD=sqrt(nv);

  

function obj=loglikb(lambda,Phi,y,nv,Z,p,type)
%LOGLIKB computes the negative loglikelihood using the kernel K=\lambdaK_B,
%        where K_B is the kernel using the basis functions
%        
%        Input: lambda >=0 scaling factor
%               Phi regression matrix
%               V=[v_1 ... v_n] matrix containing the basis functions  v_k
%               type=             'yy' BJ NP model
%                                 'ny' OE NP model 
%                                 'yn' time series
%
%        Output: obj negative loglikelihood


switch type
    case 'ny'
        B=kernel(lambda(1),lambda(2),p)*Z;
        V=B'*Phi'*Phi*B+nv*eye(size(B,2));
        yf=B'*y(:,1);
    case 'yn'
        B=kernel(lambda(1),lambda(2),p)*Z;
        V=B'*Phi'*Phi*B+nv*eye(size(B,2));
        yf=B'*y(:,1);
    case 'yy'
        B1=kernel(lambda(1),lambda(3),p)*Z(:,:,1);
        B2=kernel(lambda(2),lambda(4),p)*Z(:,:,2);
        V=[B1'*Phi(:,:,1)'*Phi(:,:,1)*B1 B1'*Phi(:,:,1)'*Phi(:,:,2)*B2; 
            B2'*Phi(:,:,2)'*Phi(:,:,1)*B1 B2'*Phi(:,:,2)'*Phi(:,:,2)*B2]+nv*eye(2*size(B1,2));
        yf=[B1'*y(:,1); B2'*y(:,2)];
end

try
    L = chol(V)';
    Li = L^-1;
    z = Li*yf; 
    obj = 0.5*(2*sum(log(diag(L))) -nv^-1*z'*z);
catch
    obj = inf;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function L= kernel(lambda,b,p)
L=sqrt(1-b)^3*sqrt(lambda)*diag(exp(-b/2*(1:p)));
