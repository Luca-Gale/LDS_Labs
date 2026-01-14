function SURE=surek(m,df)
%SUREK Computes the SURE estimator from a nonparametric model
%   SURE = SURE(Model,df)
%
%   Model = Any identfied model (IDPARAMETRIC or IDNLMODEL)
%   df = degrees of freedom of the model structure
%
%   SURE = SURE estimator V(1+ 2df/N)
%   where V is the loss function, df are the degrees of freedom of the 
%   estimated model and N is the number of estimation data.
%

%   Mattia Zorzi December 2019



p=size(m.Report.Parameters.ParVector,1);
N=m.Report.DataUsed.Length;
SURE=m.Report.Fit.LossFcn*(1+2*df/N);


 