% this resolves the unidentifiability issue with the allele frequency of
% the eQTL 

% If called to resolve a Fit from M1, P_AonR = P_AonA, and we want to know
% whether A is the High or the the Low expressed eQTL allele.

% If A turns out to be the Low expressed allele, we flip the allele
% frequency so that it represent the High expressed allele instead. 

function MLE = Resolve_Frequency_Flip_BLN(X, Xc, Xs, MLE)
%% classify based on MLE and separate the central classs
[MAP, C] = get_MAP_Class_BLN(X, Xc, MLE);
G2_Filt = C~=2;% This is those that fall under class 1 and 3 in the ASE data which are actually the heterozyguous cases in eQTL (G2) 

%% cluster the points using the known frequency and fold change
% normalize counts by Xs, add pseudo-count and log-transform!
PseudoCount = 1;
D = log(PseudoCount + Xs .* (X+Xc));

% calculate class probabilities in ASE data
p_H_on_R = MLE.Conditionals(1);
p_H_on_A = MLE.Conditionals(2);

p_HH   = p_H_on_R .* p_H_on_A;
p_HL   = p_H_on_R .* (1-p_H_on_A) + (1-p_H_on_R) .* p_H_on_A;
p_LL   = (1-p_H_on_R) .* (1-p_H_on_A);


% calculate exprected fold change between classes
s_H = Logit(MLE.P(1))-Logit(MLE.P(2));

rMu = [               % relative expression wrt A/a eQTL genotype
    log(2)+s_H        % G1: HH
    log(1+exp(s_H))   % G2: HL/LH
    log(2)            % G3: LL
    ];

P = [               % Genotype probabilities for the eQTL
    p_HH            % G1: HH
    p_HL            % G2: HL/LH
    p_LL            % G3: LL
    ];

% cluster the points using the frequency and fold change
[LL, ~, ~] = ConstrainedGaussianTriMixture(G2_Filt, D, rMu, P); 

% cluster the points using one minus the frequency and fold change
[LL_r, ~, ~] = ConstrainedGaussianTriMixture(G2_Filt, D, rMu, P(3:-1:1)); 

% report the best
if LL_r > LL
    % reverse the frequency!
    MLE.Conditionals = 1-MLE.Conditionals;
end
end
