%% Example run for a set of AE files

%% These are some basic arguments. 
% Each script has some additional parameters that are a little more involved and can be tuned from inside the code.
AnalysisLabel   = 'Test_run'; % your favorite label for this analysis
minOccurance    = 10; % aeSNPs with data in less than this many people in the entire data will be discarded right off the bat from the file "Common_SNPs".
FileNamePattern = '01_Data/Sample_AE_Data/*.txt'; % the address for all ASE file
SNP_frq_file    = '01_Data/Sample_Allelefrequency_Data/Sample_SNP_Allelefrequencies.txt'; % See inside S02_Add_AlleleFreqs.m for details.
tot_Indv_Used_THR = 5; % Minimum number of individuals available in data for a gene estimate before it's reported. For GTEx we used **50**

%% Run all scripts one after another
S01_Prepare_ASE_data(FileNamePattern, minOccurance, AnalysisLabel)
S02_Add_AlleleFreqs(AnalysisLabel, SNP_frq_file)
S03_Fit_BLNmixture(AnalysisLabel)
S04_Estimates_Vg(AnalysisLabel)
S05_Get_Genewise_Vg(AnalysisLabel, tot_Indv_Used_THR)