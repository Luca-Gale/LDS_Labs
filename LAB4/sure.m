function SURE=sure(m)
%SURE Computes the SURE estimator from a model
%   SURE = SURE(Model)
%
%   Model = Any identfied model (IDPARAMETRIC or IDNLMODEL)
%
%   SURE = SURE estimator V(1+ 2d/N)
%   where V is the loss function, d is the number of estimated parameters
%   and N is the number of estimation data.
%

%   Mattia Zorzi 2017



p=size(m.Report.Parameters.ParVector,1);
N=m.Report.DataUsed.Length;
SURE=m.Report.Fit.LossFcn*(1+2*p/N);
