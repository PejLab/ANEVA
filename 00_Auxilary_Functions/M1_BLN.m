% X is ASE ref count
% Xc is ASE alt count
% r0 is the reference ratio due to bias. 

% Optimization params:

% Z(1): Minor allele frequency of eQTL SNP assuming no LD with the ASE SNP. We can't identify if the higher expressed eQTL allele is major or the minor allele. We resolve this using total gene
% expression data in 'Resolve_Frequency_Flip_BLN.m'

% Z(2): is s_{H,L}, that is the abs. log aFC of the eQTL SNP.

% Z(3): is sigma_r, standard dev of the residual regulatory effects in log
% aFC. The std of the Normal dist. in the BLN dist.
function [MLE] = M1_BLN(X,Xc, r0, NRepeats)

s0=Logit(r0); % r0 is the reference ratio due to bias. 


OPToptions = optimoptions(@fmincon, 'MaxFunEvals', 10000, 'MaxIter', 10000, 'GradObj', 'off', 'Display', 'notify');%,  'Display', 'iter-detailed');
Problem2.objective = @(Z)BLN_mixture_Likelihood(X,Xc,...
    ThreeDensities(Z(1), Z(1)),...
    Logistic([s0+Z(2), s0, s0-Z(2)]), ...
    Z(3)*ones(1,3));
Problem2.solver = 'fmincon';
Problem2.options = OPToptions;
Problem2.UseParallel = true;
% Problem2.TolX = 1E-4;
Problem2.ScaleProblem = 'obj-and-constr';
Problem2.lb = [0 0 0];
Problem2.ub = [.5 log(100) 1.5];

for R = NRepeats:-1:1
    Problem = Problem2;
    if R==1
        Pinit = (r0+1)/2;
    else
        Pinit = rand*(1-r0)+r0;
    end
    Problem.x0 = [.25 Logit(Pinit) 1];
    
    [optZ,MLE(R).NLL,MLE(R).exitflag] = fmincon(Problem);
    
    %% Make outputs
    MLE(R).Conditionals       = ones(2,1) * optZ(1);
    MLE(R).s_H     = optZ(2);
    MLE(R).W       = ThreeDensities(MLE(R).Conditionals(1), MLE(R).Conditionals(2));
    MLE(R).P       = Logistic([s0+optZ(2), s0, s0-optZ(2)]);
    MLE(R).Std     = optZ(3)*ones(1,3);
    MLE(R).RegulatorInLD=0;
    
end

%% Choose the best initilization
[~, IminR] = min(vertcat(MLE(:).NLL));  
MLE = MLE(IminR);

%% Claculate BIC
v = 1 + 1 + 1; % P(H), s_{H,L}, and sigma_r
MLE.BIC = 2*MLE.NLL +  v * log(length(X));
end