% This script adds the allele frequencies for the aeSNPs to the data. This
% information is needed for calculating VG later on.
% Inputs:
% "AnalysisLabel" is the label given to the analysis =, should be the same as
% what was used for S01_Prepare_ASE_data.m

% "SNP_frq_file" is a text file containing two columns: snp_ids, and Ref_allele_freq
% The order of the column/rows in this file are not important. Also the file can
% have any other number of columns/rows, just the larger files will be slower to load. 
% Also the script uses the snp_ids to match rows in this file with the ASE
% data. So the snp_ids here should be the same as those used before(see ColumnLabel_for_uniqueID in S01_Prepare_ASE_data.m)


function S02_Add_AlleleFreqs(AnalysisLabel, SNP_frq_file)
if nargin == 0
    % Make a test run
    AnalysisLabel= 'Test_run';
    SNP_frq_file = ['01_Data/Sample_Allelefrequency_Data/Sample_SNP_Allelefrequencies.txt'];    
end
S01_Outputfile = ['02_Results/' AnalysisLabel '/Common_SNPs.mat'];


disp('Reading in Frqs...')
Frqs = Pej_Read_Table(SNP_frq_file);
Add_AlleleFreqs(S01_Outputfile, Frqs);
fprintf('\tdone!')

end

function Add_AlleleFreqs(File, Frqs)
try
    load(File);
   
    [~, ai, bi] = intersect(jSNPcommon.uniqueID,Frqs.snp_ids);
    
    N = length(jSNPcommon.uniqueID);
    jSNPcommon.REFfreq = nan(N,1);
    jSNPcommon.REFfreq(ai) = Frqs.Ref_allele_freq(bi);
    jSNPcommon.ALTfreq     = 1 - jSNPcommon.REFfreq;
    
    disp([num2str(length(bi)/N*100) '% matched'])
    save(File, 'jSNPcommon', 'MetaData');
    
    if (length(bi)/N)<.9
        warning on
        beep
       warning('Less than 90% of the aeSNPs were found in SNP frquency file.') 
    end
     
catch errr
    disp([File ' failed!'])
    errr.getReport
end
end
