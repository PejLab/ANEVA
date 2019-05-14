% This file gets a set of ASE data files and merges them into a matrix,
% it make a Union of all SNPs putting NaN into missing values. The folmat
% of these files is described/defined further below in the section "%% Constants/ Input file format"

% "FileNAmePattern" indicates the input files, it can be for example:

%   Myfoler/  all files in this this folder
%   Myfolder/*.txt
%   *.[1-9].dd

% Basically any sort of pattern that fits to linux find function would fit
% here too. '

% Example: S01_Prepare_ASE_data('01_Data/AE_Sample_Data/*.txt', 10, 'Test_run');

% Pej March 2015 NYGC
%--------------------
function S01_Prepare_ASE_data(FileNamePattern, minOccurance, AnalysisLabel)
if nargin == 0
    % Make a test run
    AnalysisLabel= 'Test_run';
end

OutFolder = ['02_Results/' AnalysisLabel];
mkdir(OutFolder);

delete([OutFolder '/S01_Prepare_ASE_data_run.log']);
diary( [OutFolder '/S01_Prepare_ASE_data_run.log']);

disp('Run started at:')
disp(datetime(now,'ConvertFrom','datenum'))

if nargin == 0
    beep
    disp('No input parameters were given. Running demo on sample data.')
    minOccurance = 10;
    FileNamePattern = '01_Data/AE_Sample_Data/*.txt';
end

%% Run the script
try
    S01_Prepare_ASE_data_(FileNamePattern, minOccurance, OutFolder)
catch ee
    ee.getReport
    diary off
end

end

function S01_Prepare_ASE_data_(FileNamePattern, minOccurance, OutFolder)
%% Constants/ Input file format (PARAMETERS)

% Mandatory fields ======================
ColumnLabel_for_uniqueID   = 'VARIANT_ID';  % Label for the column carrying a unique varinat ID for aeSNP.      **Edit this to match your data**
ColumnLabel_for_refCount   = 'REF_COUNT';   % Label for the column carrying ref. allele count at aeSNP.         **Edit this to match your data**
ColumnLabel_for_altCount   = 'ALT_COUNT';   % Label for the column carrying alt. allele count at aeSNP.         **Edit this to match your data**
ColumnLabel_for_GeneID     = 'GENE_ID';     % Label for the column carrying the gene ID associated with aeSNP.  **Edit this to match your data**


% Optional fields (if you don't need these comment out the lines)
ColumnLabel_for_null_ratio = 'NULL_RATIO';  % Label for the column carrying the null ratio for the aeSNP. This is the reference bias, we define it as the expected reference ratio in absence of any regulatory difference between the haplotype. If you don't have this in your data, delete this line and we assume null ratio is 0.50  **Edit this to match your data**

% Optional QC columns (if you don't need these comment out the lines)
% if there are columns "MAPPING_BIAS_SIM", "GENOTYPE_WARNING", and
% "LOW_MAPABILITY" are available in the data, all lines with 1 in either of
% those will be removed before the analysis too. If you have some other QC
% things just apply to the data before running ANEVA.
QC_ColumnLabels = {'MAPPING_BIAS_SIM', 'GENOTYPE_WARNING', 'LOW_MAPABILITY'}; % QC_ColumnLabels is a cell array of column names for qulity control columns in data.


MinCount = 8; % Discard ASE data with total count smaller than this from the data
%% Default values
Files = Pej_GetFiles(FileNamePattern);
NF = length(Files);
if nargin <2
    minOccurance = min(10, NF); % Variants that are less common than this will be discarded from the file "Common_SNPs".
end

%% The rest of the code
fprintf(1, '%d files found!\n', NF);
for i = NF:-1:1
    fprintf(1, '%i files left...', i);
    
    %     FormatS = '%s%s%s%s%s%s%s%s%f%f%f%f%f%f%f%f%s%s%s%s%s%f%f%f'; % for standard GTEx_V6
    %     FormatS='%s%s%s%s%s%f%f%f%f%f%f%f%f%s%s%f%s%s%f%f%f%f%f%f'; % For GEUVADIS
    
    %     tmpD= Pej_Read_Table(Files{i}, [], false, FormatS);
    tmpD= Pej_Read_Table(Files{i}, [], false, [], true);
    % converting the data column names to what is used around the code.
    D{i}.uniqueID   =tmpD.(ColumnLabel_for_uniqueID);
    D{i}.refCount   = str2double(tmpD.(ColumnLabel_for_refCount));
    D{i}.totalCount = D{i}.refCount + str2double(tmpD.(ColumnLabel_for_altCount));
    D{i}.GeneID     =tmpD.(ColumnLabel_for_GeneID);
    %     D{i}.chr        = tmpD.CHR;
    %     D{i}.position   = tmpD.POS;
    %     D{i}.Variant_Annotation     = repmat({'NA'}, size(D{i}.GeneID,1), size(D{i}.GeneID,2));
    
    
    MetaData.SAMPLE_ID{i}   = tmpD.SAMPLE_ID{1};
    %     MetaData.SUBJECT_ID{i}  = tmpD.SAMPLE_ID{1}(1:find(tmpD.SAMPLE_ID{1}=='.',1)-1);
    
    %% If null ratio is missing from the data, leave it empty
    try
        D{i}.null_ratio = str2double(tmpD.(ColumnLabel_for_null_ratio));
    catch
        D{i}.null_ratio = nan(size(D{i}.uniqueID));
    end
    
    %% Discard low quality aeSNPs if the columns are there, if they are not it'll just skip it.
    % See Castel et al. 2015 for reference (https://doi.org/10.1186/s13059-015-0762-6)
    tmpFilt = false(size(D{i}.uniqueID));
    for i_QC = 1:length(QC_ColumnLabels)
        try
            tmpFilt = tmpFilt | str2double(tmpD.(QC_ColumnLabels{i_QC}))==1;
        catch err
            error(['Optional QC column "' QC_ColumnLabels{i_QC} '" is not available in the data. If it is supposed to be there check if there is a typo. If you don''t have or need it remove the name of this column from the list of QC columns at the begining of S01_Prepare_ASE_data.m'])
        end
    end
    fprintf('\n%d our of %d (%.2f%%) variants were filtered from %s\n', sum(tmpFilt), length(tmpFilt), (sum(tmpFilt)/length(tmpFilt))*100, tmpD.SAMPLE_ID{1});
    D{i}= Pej_Struct_RowDel(D{i}, tmpFilt);
    
    D{i} = Pej_Struct_RowSelect(D{i}, D{i}.totalCount>=MinCount);
end


mkdir(OutFolder);
% save([OutFolder '/All_Imported'], 'D', 'FileNamePattern');
% load('All_Imported');

%% Join them in a single matrix
jSNPs = Union_stats(D);
figure
hist(jSNPs.occurance(jSNPs.occurance>1), 1:max(jSNPs.occurance));
% xlim([0 101])
xlabel('Hetrozygous occurance in dataset')
ylabel('Number of aeSNPs')
Pej_SavePlot(gcf, [OutFolder '/Figures/aeSNP_occurrence']);

%%
jSNPcommon= Pej_Struct_RowDel(jSNPs, jSNPs.occurance<minOccurance); clear jSNPs
jSNPcommon.totalCount = nan(length(jSNPcommon.uniqueID),NF);
jSNPcommon.refCount   = nan(length(jSNPcommon.uniqueID),NF);
jSNPcommon.null_ratio = nan(size(jSNPcommon.uniqueID));
jSNPcommon.GeneID     = cell(size(jSNPcommon.uniqueID));
% jSNPcommon.chr        = cell(size(jSNPcommon.uniqueID));
% jSNPcommon.position   = cell(size(jSNPcommon.uniqueID));
% jSNPcommon.Variant_Annotation = cell(size(jSNPcommon.uniqueID));

for i = 1:NF
    [~, Ij, Id] = intersect(jSNPcommon.uniqueID, D{i}.uniqueID);
    jSNPcommon.totalCount(Ij,i) = D{i}.totalCount(Id);
    jSNPcommon.refCount(Ij,i)   = D{i}.refCount(Id);
    %     jSNPcommon.chr(Ij,1)        = D{i}.chr(Id); % This line does not really have o be done so many times but it's a lazy apporaoch to make sure all SNPs are covered.
    %     jSNPcommon.position(Ij,1)   = D{i}.position(Id); % This line does not really have o be done so many times but it's a lazy apporaoch to make sure all SNPs are covered.
    jSNPcommon.null_ratio(Ij,1) = D{i}.null_ratio(Id); % This line does not really have o be done so many times but it's a lazy apporaoch to make sure all SNPs are covered.
    jSNPcommon.GeneID(Ij,1)     = D{i}.GeneID(Id); % This line does not really have o be done so many times but it's a lazy apporaoch to make sure all SNPs are covered.
    %     jSNPcommon.Variant_Annotation(Ij,1)     = D{i}.Variant_Annotation(Id); % This line does not really have o be done so many times but it's a lazy apporaoch to make sure all SNPs are covered.
    
    % delete it
    D{i} = [];
end


%% Read in sample expression normalization factors
save([ OutFolder '/Common_SNPs'], 'jSNPcommon', 'MetaData');%,'-v7.3'); % ver 7.3 allows for saving large arrays.
disp(['Data for ' num2str(length(jSNPcommon.uniqueID))    ' aeSNPs was saved in: '  OutFolder '/Common_SNPs.mat'] )

tfout = fopen([OutFolder '/Common_SNPs_SNPids.txt'], 'w');
fprintf(tfout, 'snp_ids\n');
fprintf(tfout, '%s\n', jSNPcommon.uniqueID{:});
fclose(tfout);
disp(['Used aeSNP IDs were saved in: ' [OutFolder '/Common_SNPs_SNPids.txt']] )


disp('Run ended at:')
disp(datetime(now,'ConvertFrom','datenum'))
diary off
end


function U = Union_stats(D)
tmpFile = ['temp_' num2str(rand*1E+16) '_' num2str(cputime)];
tmpF = fopen(tmpFile, 'w');
for i = 1:length(D)
    fprintf(tmpF, '%s\n', D{i}.uniqueID{:});
end
fclose(tmpF);
system(sprintf('sort %s >%s.sorted', tmpFile, tmpFile));
system(sprintf('uniq -c %s.sorted >%s.counts', tmpFile, tmpFile));
system(sprintf('rm %s', tmpFile));
system(sprintf('rm %s.sorted', tmpFile));

tmpF = fopen(sprintf('%s.counts', tmpFile), 'r');
tmpIDs = textscan(tmpF, '%f%s');
fclose(tmpF);
U.uniqueID  = tmpIDs{2};
U.occurance = tmpIDs{1};
system(sprintf('rm %s.counts',tmpFile));
end


