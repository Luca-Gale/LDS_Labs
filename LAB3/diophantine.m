function [D, N] = diophantine(A, B, D0)
% DIOPHANTINE   Get a solution of the diophantine equation.
%
%   [D,N] = DIOPHANTINE(A,B,D0) get the min deg solution 
%   wrt N(z) of the following diophantine eqn:
%
%       D(z)A(z) + N(z)B(z) = D0(z)
%
%
 

%%  Define the degs of D and N (ie k and h)
n = length(A)-1;
m = length(B)-1;
l = length(D0)-1;

if l >= m+n,
    k = l-n;
else
    k = m-1;
end;
h = n-1;

%%  Define the Sylvester matrix associated with A(z) and B(z)
M = zeros(k+1+n, k+1+n);
for i=1:k+1,
    for j=0:n,
        M(i+j,i) = A(n+1-j);
    end;
end;

for i=1:n,
    for j=0:m,
        M(i+j,i+k+1) = B(m+1-j);
    end;  
end;

%%  Solve the Diophantine equation
C = zeros(k+1+n,1);
for i=1:l+1,
    C(i) = real(D0(l+2-i));
end;
X = pinv(M)*C;

D = flipud(X(1:k+1)).';
N = flipud(X(k+2:k+2+h)).';

end
