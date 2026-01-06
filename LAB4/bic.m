function BIC1 = bic(model)
%BIC Computes B-Akaike's Information Criterion(BIC) from a model
%   BIC = BIC(Model)
%
%   Model = Any identfied model (IDPARAMETRIC)
%
%   BIC = B-Akaikes Information Criterion V(1 + dlogN/N)
%   where V is the loss function, d is the number of estimated parameters
%   and N is the number of estimation data.
%
%

%  Mattia Zorzi, December 2019


V = model.EstimationInfo.LossFcn;
N = model.EstimationInfo.DataLength;
d = size(model.report.parameters.ParVector,1);

BIC1 =log(V)+d*log(N)/N;
