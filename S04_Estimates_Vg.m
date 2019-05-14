% Example: S04_Estimates_Vg({'01_Data/Synthetic_ASE_Data/Common_SNPs.mat'})

function S04_Estimates_Vg(AnalysisLabel)
addpath(['00_Auxilary_Functions/'])

if nargin ==0
    % Make a test run
    AnalysisLabel= 'Test_run';
end
DataPath = ['02_Results/' AnalysisLabel '/Common_SNPs.mat'];

% Precalculate the function
PreCalculated =  './00_Auxilary_Functions/S04_Calculate_stats_Lookup_MAD_BLN_v3.mat';
if exist(PreCalculated, 'file')
    warning('%s exists, integral will be imputed from this file, if you want it to be fresh calculated DELETE the file and rerun!', PreCalculated)
else
    Make_Dg_Stats_Lookup(PreCalculated)
end

Lookup = load(PreCalculated);
get_Dg_Estimates(DataPath, Lookup)

end

function get_Dg_Estimates(DataPath, Lookup)
% =========================
% PARAMETERS
Mincases=6; % Only analyze genes with "Mincases" times have expression above "MinExpression"
MinExpression=30; % Only analyze genes with expression higher than "MinExpression" in at least "(Mincases)" samples
% =========================
% Pejman Feb 2016
% =========================
addpath(['00_Auxilary_Functions/'])

[DataFolder,FileName,~] = fileparts(DataPath);
N0=0;KN=1;
ResultFile = [DataFolder '/' FileName '_stats_BLN_v3'];
S03_Fit_Model_OutDump    = [DataFolder '/' FileName '_Constrained_3Clusters_BLN.mat'];

load(S03_Fit_Model_OutDump); % get the MLE results
load(DataPath); % get original data
snp = jSNPcommon; clear jSNPcommon jSNPs

N = length(snp.uniqueID);
usedFilt = false(N,1);
snp.Dg = nan(N,1);
snp.Vg = nan(N,1);

snp.Descriptions = cell(N,1);

M0_F = false(N,1);

tic
for i = (N-N0):-KN:1
    if i/5000 == round(i/5000)
        fprintf('%d ', i);
    end
    
    try
        nr = snp.null_ratio(i);
        if isnan(nr); nr = .5;end
        
        F = ~isnan(snp.totalCount(i,:));
        tN = snp.totalCount(i,F);
        
        if sum(tN>=MinExpression)>=Mincases && ~isempty(MLE(i,1).NLL)
            usedFilt(i)= true;
            
            S = Logit(MLE(i).P(1))-Logit(MLE(i).P(2));
            if S==0
                %                 it's M0
                M0_F(i) = true;
            else
                p_H = MLE(i).Conditionals(1).*snp.REFfreq(i) + MLE(i).Conditionals(2).*snp.ALTfreq(i); % pop frequency of the Higher expressed eQTL allele
                [snp.Dg(i), snp.Vg(i)] = Expected_Expression_MAD_BLN_v3(S, MLE(i).Std(1), p_H);
            end
        end
    catch EEEE
        %         discard
        snp.Descriptions{i} = sprintf('SNP%d\t%s\tFailed: %s\n', i, snp.uniqueID{i}, EEEE.message);
        EEEE.getReport
    end
end

tODs = vertcat(MLE(M0_F).Std);
snp.Dg(M0_F)  = interp1(Lookup.sigma_r, Lookup.Dg, tODs(:,1));
snp.Vg(M0_F)  = interp1(Lookup.sigma_r, Lookup.Vg, tODs(:,1));

OutBuff.uniqueID            = snp.uniqueID(usedFilt);
OutBuff.GeneID              = snp.GeneID(usedFilt);
OutBuff.Dg                  = snp.Dg(usedFilt);
OutBuff.Vg                  = snp.Vg(usedFilt);
OutBuff.Descriptions        = snp.Descriptions(usedFilt);
try
    OutBuff.Variant_Annotation  = snp.Variant_Annotation(usedFilt);
catch
end


snp = Pej_Struct_RowSelect(snp, usedFilt);
save(ResultFile, 'snp');
Pej_Write_Table([ResultFile '.txt'], OutBuff);
toc
end


