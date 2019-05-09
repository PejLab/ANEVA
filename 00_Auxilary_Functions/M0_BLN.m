% X is ASE ref count
% Xc is ASE alt count
% r0 is the reference ratio due to bias. 

% optimization params:
% Z(1): is sigma_r, standard dev of the residual regulatory effects in log
% aFC. The std of the Normal dist. in the BLN dist.

function MLE = M0_BLN(X,Xc, r0)
NRepeats = 1;
s0=Logit(r0); % r0 is the reference ratio due to bias. 

OPToptions = optimoptions(@fmincon, 'MaxFunEvals', 10000, 'MaxIter', 10000, 'GradObj', 'off', 'Display', 'notify');%,  'Display', 'iter-detailed');
Problem2.objective = @(Z)BLN_mixture_Likelihood(X,Xc,...
    1,...
    r0, ...
    Z(1));

Problem2.solver = 'fmincon';
Problem2.options = OPToptions;
Problem2.UseParallel = true;
% Problem2.TolX = 1E-4;
Problem2.ScaleProblem = 'obj-and-constr';
Problem2.lb = 0;
Problem2.ub = 1.5; % 

% parfor R = 1:NRepeats
for R = NRepeats:-1:1
    Problem = Problem2;
    Problem.x0 = 1;
    
    [optZ,MLE(R).NLL,MLE(R).exitflag] = fmincon(Problem);
    
    %% Make outputs
    MLE(R).Conditionals       = zeros(1,2);
    MLE(R).s_H     = 0;
    MLE(R).W       = ThreeDensities(MLE(R).Conditionals(1), MLE(R).Conditionals(2));
    MLE(R).P       = Logistic([s0, s0  ,s0]);
    MLE(R).Std     = optZ(1)*ones(1,3);
    MLE(R).RegulatorInLD= -1; % The null model with only residual variation
    
end

%% Choose the best initilization
[~, IminR] = min(vertcat(MLE(:).NLL));
MLE = MLE(IminR);

%% Claculate BIC
v = 1; % std
MLE.BIC = 2*MLE.NLL +  v * log(length(X));
end

