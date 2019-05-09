% This function calculates Vg and Dg caused by a single eQTL using
% simulation for one signle input set


% Dg: Mean Abs Dev in *natural* log scale due to genetic effects
% Vg: Variance     in *natural* log scale due to genetic effects

function [Dg, Vg] = Simulate_Expected_Expression_MAD_byeQTL(log2_aFC_A2R, fR)
N = 10000; % individuals
H1geno = rand(N,1)>fR; % 1means aleternative
H2geno = rand(N,1)>fR; % 1means aleternative

aFC_A2R = 2.^log2_aFC_A2R;
E0 = 1; % basal exp in abosolute scale

H1exp = E0 * ones(N,1);
H1exp(H1geno==1) = H1exp(H1geno==1) * aFC_A2R;

H2exp = E0 * ones(N,1);
H2exp(H2geno==1) = H2exp(H2geno==1) * aFC_A2R;

TotExp = H1exp + H2exp ;

logTotExp = log(TotExp);
Vg = var(logTotExp);
Dg = mad(logTotExp);
end