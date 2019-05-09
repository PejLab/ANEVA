%  This classifies AE data into classes based on a mixture of
%  beta-binomials parameters given in MLE

function [MAP, C] = get_MAP_Class_BLN(X, Xc, MLE)
K = length(MLE.P); % number of classes
LW = (MLE.W /sum(MLE.W ));
for i = K:-1:1
    Mu = Logit(MLE.P(i))-Logit(MLE.P(2));
    L(:,i) = Pej_pdf_BLN(X, Xc, Mu, MLE.Std(i).^2).*LW(i);
end
L = L./repmat(sum(L,2),1,K); % posterior
[MAP, C] = max(L,[],2);

end