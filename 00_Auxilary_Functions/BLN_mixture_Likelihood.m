function NLL = BLN_mixture_Likelihood(X,Xc,W,P,Std)
% Calculate negative total log-likelihood of a Binomial Logit Normal mixture
% X is an array of counts corresponding to positives cases
% Xc is an array of counts corresponding to negative cases

% W is the mixture wieghts, it should sum up to 1.
% P 0s an array of fractions for the Normal expectation
% STD is an array controling the variance
% Pej March 2015
%----------
LW = W/sum(W);
Mu = Logit(P);
V  = Std.^2;
for i = length(W):-1:1
    L(:,i) = Pej_pdf_BLN(X, Xc, Mu(i), V(i), false)*LW(i);
end
NLL = -sum(log(sum(L,2)));%+abs(sum(W)-1);
end
