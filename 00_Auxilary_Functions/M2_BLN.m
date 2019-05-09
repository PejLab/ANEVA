% X is ASE ref count
% Xc is ASE alt count
% r0 is the reference ratio due to bias. 

% Optimization params:

% Z(1): P(H|R), that is probability of observing higher expressed eQTL SNP
% on the same hap. as the Ref. ASE allele.

% Z(2): P(H|A), that is probability of observing higher expressed eQTL SNP
% on the same hap. as the Alt. ASE allele.

% Z(3): is s_{H,L}, that is the abs. log aFC of the eQTL SNP.

% Z(4): is sigma_r, standard dev of the residual regulatory effects in log
% aFC. The std of the Normal dist. in the BLN dist.

function MLE = M2_BLN(X,Xc, r0, NRepeats, InitParams)

s0=Logit(r0); % r0 is the reference ratio due to bias. 

OPToptions = optimoptions(@fmincon, 'MaxFunEvals', 10000, 'MaxIter', 10000, 'GradObj', 'off', 'Display', 'notify');%,  'Display', 'iter-detailed');
Problem2.objective = @(Z)BLN_mixture_Likelihood(X,Xc,...
    ThreeDensities(Z(1), Z(2)),...
    Logistic([s0+Z(3), s0, s0-Z(3)]), ...
    Z(4)*ones(1,3));
Problem2.solver = 'fmincon';
Problem2.options = OPToptions;
Problem2.UseParallel = true;
% Problem2.TolX = 1E-4;
Problem2.ScaleProblem = 'obj-and-constr';
Problem2.lb = [0 0 0 0];%[0 0  logitZ2 log(2)   ];
Problem2.ub = [1 1 log(100)   1.5]; % 8 is as far as usually you go from technical repeats

for R = NRepeats:-1:1
    Problem = Problem2;
    if R==1
        Pinit = (r0+1)/2;
    else
        Pinit = rand*(1-r0)+r0;
    end
    Problem.x0 = [.5 .5 Logit(Pinit) 1];
    
    if R==1 && nargin>4 && ~isempty(InitParams)
        %% Override the initial parameters
        Problem.x0 = InitParams;
    end
    
    
    [optZ,MLE(R).NLL,MLE(R).exitflag] = fmincon(Problem);
    
    %% Make outputs
    MLE(R).Conditionals       = optZ(1:2);
    MLE(R).s_H     = optZ(3);
    MLE(R).W       = ThreeDensities(MLE(R).Conditionals(1), MLE(R).Conditionals(2));
    MLE(R).P       = Logistic([s0+optZ(3), s0, s0-optZ(3)]);
    MLE(R).Std     = optZ(4)*ones(1,3);
    MLE(R).RegulatorInLD=1;
end

%% Choose the best initilization
[~, IminR] = min(vertcat(MLE(:).NLL));
MLE = MLE(IminR);

%% Claculate BIC
v = 2 + 1 + 1; % P(H|R), P(H|A), s_{H,L}, and sigma_r
MLE.BIC = 2*MLE.NLL +  v * log(length(X));
end