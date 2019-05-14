% This file collapses the SNP level estimates to produce gene level ones,
% and also makes a file that includes all genewise estimates 

% tot_Indv_Used_THR % Minimum number of individuals available in data for a gene estimate before it's reported. For GTEx we used 50

function S05_Get_Genewise_Vg(AnalysisLabel, tot_Indv_Used_THR)
%=========================
% PARAMETERS
MinExpression   = 30; % Only consider cases were at least one aeSNP is covered by more than this number of reads in an individual for tot_Indv_Used_THR
%=========================

if nargin ==0
    % Make a test run
    AnalysisLabel= 'Test_run';
    tot_Indv_Used_THR = 5;
end
DataPath = ['02_Results/' AnalysisLabel '/Common_SNPs_stats_BLN_v3.mat'];

load(DataPath);
RT=snp; clear snp
%% Handle multiple SNPs for one gene
% Average over all SNPs
[~,I] = sort(RT.GeneID);
RT = Pej_Struct_RowSelect(RT,I);

SZ = size(RT.uniqueID);
CurrGene = RT.GeneID(1); i0=1; iu=0;
isUniq=false(SZ);isUniq(1)=true;
RT.Dg_GeneWise = nan(SZ);
RT.Vg_GeneWise = nan(SZ);
RT.tot_Indv_Used = nan(SZ);

for i = 2:length(RT.uniqueID)
    if strcmp(CurrGene, RT.GeneID(i)) && ~strcmp('NaN', CurrGene)
        % same gene
        
    else
        % new gene or un-annotated gene
        isUniq(i) = true;
        iu = iu+1; % number of unique hits
        
        gi=(i0:i-1);
        RT.tot_Indv_Used(gi) = sum(any((RT.totalCount(gi,:)>=MinExpression) & repmat(~isnan(RT.REFfreq(gi,:)), 1,size(RT.totalCount,2)),1));
        
        tmpM = Pej_Mean_withNaNs(RT.totalCount(gi,:),2).*sum(RT.totalCount(gi,:)>0,2); % total number of reads that contributed to this estimate
        RT.Dg_GeneWise(gi)      = wSummary(RT.Dg(gi), tmpM);
        RT.Vg_GeneWise(gi)      = wSummary(RT.Vg(gi), tmpM);
        % update current gene
        CurrGene = RT.GeneID(i);
        i0 = i;
    end
end

RT = Pej_Struct_RowSelect(RT, isUniq);

RT_OUT.GeneID = RT.GeneID;
RT_OUT.Dg_GeneWise = RT.Dg_GeneWise;
RT_OUT.Vg_GeneWise = RT.Vg_GeneWise;

RT_OUT = Pej_Struct_RowSelect(RT_OUT, RT.tot_Indv_Used>=tot_Indv_Used_THR);

[DataFolder,FileName,~] = fileparts(DataPath);
Pej_Write_Table([DataFolder '/' FileName '-GeneWise.txt' ], RT_OUT);

end

function Y = wSummary(X,W)
t = X;
tw = W;
F = ~isnan(t+tw);
if sum(F)>=1  % THIS IS A TUNABLE PARAMETER! THE MINIMUM NUMBER OF SCORES
    Y = 1./(sum((1./t(F)).*tw(F))/sum(tw(F))); % Weighted harmonic mean
else
    Y=nan;
end
end