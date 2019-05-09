function [MLE] = Fit_BLN(X,Xc, null_ratio, ModelType, NRepeats, InitParams)
if nargin < 6
    InitParams = [];
end

switch ModelType
    case 2
        MLE = M2_BLN(X,Xc, null_ratio, NRepeats, InitParams);
    case 1
        MLE = M1_BLN(X,Xc, null_ratio, NRepeats);
    case 0
        MLE = M0_BLN(X,Xc, null_ratio);
end

end