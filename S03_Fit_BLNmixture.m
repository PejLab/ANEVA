function S03_Fit_BLNmixture(AnalysisLabel, N0, KN)
%=========================
% PARAMETERS
Mincases        = 6;%6; % Only analyze genes with "Mincases" times have expression above "MinExpression"
MinExpression   = 30; % Only analyze genes with expression higher than "MinExpression" in at least "(Mincases)" samples
NRepeats        = 1; % Number of times the clustering fit is done before the best one is reported. Each round get a new random initialization.
LRT_p_Thr       = 0.01; % Threshold for Likelihood ratio Test model selection

%=========================
% Pejman Sep 2017
%=========================
if nargin ==0
    % Make a test run
    AnalysisLabel= 'Test_run';
end

addpath(['00_Auxilary_Functions/'])
DataPath = ['02_Results/' AnalysisLabel '/Common_SNPs.mat'];

[DataFolder,FileName,~] = fileparts(DataPath);

if nargin<3
    N0=0;KN=1;
    OutDump    = [DataFolder '/' FileName '_Constrained_3Clusters_BLN.mat'];
else
    OutDump    = sprintf('%s/%s_Constrained_3Clusters_BLN_%dK_%d.mat', DataFolder, FileName, KN, N0);
end

if exist(OutDump, 'file')
    % This has been analyzed before
    warning('Result file, %s, was already found!\nI will skip those case that already have an MLE solution and will only analyze those that have no results in the file already. \nIf you want a complete fresh run delete the file and rerun this script.', OutDump);
    load(OutDump);
end


load(DataPath);
snp = jSNPcommon; clear jSNPcommon jSNPs
%------------------------
% Learn library normalization from ASE data. If you have better estimates
% for this from somewhere else, like from full transcriptome counts, you
% can substitute it here. 
tLibSz = Pej_Mean_withNaNs(snp.totalCount);
MetaData.Normfactor = mean(tLibSz)./tLibSz;
% -----------------------

snp.altCount=snp.totalCount-snp.refCount; snp = rmfield(snp,'totalCount');

N = length(snp.uniqueID);
ts = 0;
tic
for i = (N-N0):-KN:1
    try
        F = ~isnan(snp.altCount(i,:));
        X  = snp.refCount(i,F)';
        Xc = snp.altCount(i,F)';
        Xs = MetaData.Normfactor(F);
        nr = snp.null_ratio(i);
        if isnan(nr)
            nr = .5;
        end
        
        tN = X+Xc;
        NThrSamples = sum(tN>=MinExpression);
        if ~isnan(snp.REFfreq(i)) && NThrSamples>=Mincases && ...  % we have the allele frequency for aeSNP, and it has enough samples AND
                (~exist('MLE', 'var') || i>length(MLE) || isempty(MLE(i).BIC) || isnan(MLE(i).BIC))  % it is not already available in the results!
            disp(i)
            tMLE(2) = Fit_BLN(X, Xc, nr, 1, NRepeats);
            tMLE(2) = Resolve_Frequency_Flip_BLN(X', Xc', Xs, tMLE(2));
            
            tMLE(1) = Fit_BLN(X, Xc, nr, 0, NRepeats);
            
            InitParams = [tMLE(2).Conditionals' Logit(tMLE(2).P(1))-Logit(tMLE(2).P(2)) tMLE(2).Std(1)];
            tMLE(3) = Fit_BLN(X, Xc, nr, 2, 1, InitParams);
            
            % BIC Model selction
            %             [~, bestI] = min(vertcat(tMLE(:).BIC));
            
            % LRT Model selection
            LRT_p_M1vsM0 = chi2cdf(2 * (tMLE(1).NLL - tMLE(2).NLL), 2, 'upper'); % comparing M1 to M0
            LRT_p_M2vsM0 = chi2cdf(2 * (tMLE(1).NLL - tMLE(3).NLL), 3, 'upper'); % comparing M2 to M0
            LRT_p_M2vsM1 = chi2cdf(2 * (tMLE(2).NLL - tMLE(3).NLL), 1, 'upper'); % comparing M2 to M1
            
            if min(LRT_p_M2vsM0, LRT_p_M1vsM0) > LRT_p_Thr
                % M0 cannot be rejected
                bestI = 1;
            elseif max(LRT_p_M2vsM0, LRT_p_M1vsM0) <= LRT_p_Thr
                % Both M2 and M1 are better than M0
                if LRT_p_M2vsM1 <= LRT_p_Thr
                    % choose M2
                    bestI = 3;
                else
                    % choose M1
                    bestI = 2;
                end
            else
                %Only one of them is better than M0, pick that one
                if LRT_p_M2vsM0 <= LRT_p_Thr
                    bestI = 3;
                else
                    bestI = 2;
                end
            end
            
            bMLE = tMLE(bestI);
            bMLE.NThrSamples = NThrSamples;
            MLE(i,1) =  bMLE;
            
            ts = ts+1;
            if ts/250 == round(ts/250)
                isIntermediateOutput = true;
                save(OutDump, 'MLE', 'isIntermediateOutput');
            end
            
        else
            if exist('MLE', 'var')==0 || i>length(MLE) || isempty(MLE(i,1).NLL)
                MLE(i,1)=Fillnan();
            end
        end
    catch EEEE
        EEEE.getReport
        warning('Error Occured!')
    end
end
toc
if length(MLE)<length(snp.refCount)
    MLE(length(snp.refCount)).exitflag=nan;
end

isIntermediateOutput = false;
save(OutDump, 'MLE', 'isIntermediateOutput');
end


function MLE = Fillnan()
MLE.NLL = nan;
MLE.exitflag = nan;
MLE.Conditionals = nan(1,3);
MLE.s_H = nan;
MLE.W = nan(1,3);
MLE.P = nan(1,3);
MLE.Std = nan;
MLE.RegulatorInLD = nan;
MLE.BIC = nan;
MLE.NThrSamples = nan;
end

