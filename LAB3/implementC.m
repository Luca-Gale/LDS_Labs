function out =implementC(in)
%IMPLEMENTC implements the designed compensator


%% Get data
ep=in(1:2); % past values of the error
up=in(3:4); % past values of u
e=in(5);    % error at the present 
N=in(6:8);  % numerator of C(z)
D=in(9:11); % denominator of C(z)


%% Anti-windup paramters
uH=3;
uL=-3;

%%  Get control (include anti-windup)
v= (-D(2:end)'*up +N'*[e; ep])/D(1);
u = sat(v, uH, uL);


%%  Save "previous" vars (controller state)
ep = [e; ep(1:end-1)];
up = [u; up(1:end-1)];


%% Output
out=[ep; up];

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function usat = sat(u, uH, uL)
    if u < uL,
        usat = uL;
    elseif u > uH,
        usat = uH;
    else
        usat = u;
    end;
end