% This function gives the loglikelihood of the best Gaussian-mixture model
% with 3 components!
% fit to data, coonstrained by relative difference of the means rM3(3x1),
% and the mixing probability vector P(3x1)

% The three classes are: AA, aA/Aa, and aa

% input: 
% C2: is a binary vector marking those samples that belong to the second
% class. These pooint will not be clustered, but will affect the estimates
% and likelihood. Also No other points can be part of the second cluster.
 
% D: is a verctor of log normalized expressions

% OUtput:
% Std: is Standard deviation
function [LL, M, Std] = ConstrainedGaussianTriMixture(C2, D, rM3, P)
if ~isvector(D)
    error('INput data should be a 1xN vector, it is not a vector!')
end

if size(D,2)==1
    D = D';
end

MaxRounds = 1000; % maxumum optimization rounds


% Calculate Means 
Mu    = mean(D); % obseved mean
rM0 = sum(rM3.*P); % this is where I expect to see the mean
M = rM3 - (rM0 - Mu); % initial estimate of Mu

% Caculate Standard Dev.
Round = 0;
Std = std(D); % initial standard Dev of components

noMaxIt = true;
notCoverged = true; % relative change in data loglikelihood
RelTol = 1E-3;

Mold = M;%[nan nan];
Stdold = Std;


try
while notCoverged && noMaxIt
    Round = Round+1;
    C = Classify(C2, D, M, Std, P);
    
    
    classMeans(1,:)  = M(C);  
    M = M - mean(classMeans-D);
    
    classMeans(1,:)  = M(C);
    Std = std(classMeans-D);
    
    notCoverged = ...
        abs(Std-Stdold)/Std > RelTol || ...
        norm((M-Mold)./M) > RelTol;
    
    Mold = M;
    Stdold = Std;
    
    noMaxIt = Round<MaxRounds;
end
% disp([num2str(Round) ' rounds of optimization.'])
LL = loglike_GM(C2, D, M, Std, P);
catch E
    E.getReport
end
end

function LL = loglike_GM(C2, D, M, Std, P)
C1_3 = ~C2;
P_1 = zeros(size(D));
P_2 = zeros(size(D));
P_3 = zeros(size(D));

P_1(C1_3) = (normpdf(D(C1_3), M(1), Std)) * (P(1));
P_2(C2  ) = (normpdf(D(C2  ), M(2), Std)) * (P(2)); % those fixed to be in C2 do not mix with others
P_3(C1_3) = (normpdf(D(C1_3), M(3), Std)) * (P(3));

LL = sum(log(P_1 + P_2 + P_3));
end

function C = Classify(C2, D, M, Std, P)
C1_3 = ~C2;
P_1 = zeros(size(D));
P_3 = zeros(size(D));

P_1(C1_3) = log(normpdf(D(C1_3), M(1), Std)) + log(P(1));
P_3(C1_3) = log(normpdf(D(C1_3), M(3), Std)) + log(P(3));

C               = ones(size(D));
C(P_1 < P_3)    = 3;
C(C2)           = 2;
end