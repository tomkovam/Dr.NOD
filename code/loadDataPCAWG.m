function [tableSamplesPCAWG, tableDonorsPCAWG, nSamplesPCAWG, tableExpressionGenes, matExpressionGenesSamplesPCAWG, sSignatures, tableDriverMutationsPCAWG] = loadDataPCAWG(runAgain, expressionType, sProperties)
% Here we load all PCAWG WGS samples, link them to their RNA aliquot IDs, load expression and mutational signatures data

if (~exist('runAgain', 'var'))
    runAgain = false;
end
if (~exist('expressionType', 'var'))
    expressionType = 'fpkm_uq';
end

if (~exist('sProperties', 'var'))
    INFILE1 = 'data//PCAWG/pcawg_sample_sheet.txt';
    INFILE2 = ['data/PCAWG/expression/joint_',expressionType,'.txt'];
    INFILE3 = 'data/genes/sortedGenes.allGenes.hg19v19.txt';
    INDIR4 = 'data/PCAWG/signatures/sigProfiler/'; % INFILE4 = [DIR_DATA_ORIG, 'PCAWG/signatures/sigProfiler/PCAWG_sigProfiler_',signatureType,'_signatures_in_samples.csv'];
    INFILE5 = 'data/genes/TableS3_panorama_driver_mutations_ICGC_samples.controlled.txt';
else
    INFILE1 = sProperties.PCAWG_SAMPLE_SHEET;                                       % data/PCAWG/pcawg_sample_sheet.txt
    INFILE2 = [sProperties.PCAWG_EXPRESSION_DIR, 'joint_',expressionType,'.txt'];   % data/PCAWG/expression/joint_fpkm_uq.txt
    INFILE3 = sProperties.GENES_GENCODE;                                            % [DIR_DATA_ORIG, 'genes/sortedGenes.allGenes.hg19v19.txt'];
    INDIR4 = sProperties.PCAWG_SIGNATURES_DIR;                                      %[DIR_DATA_ORIG, 'PCAWG/signatures/sigProfiler/']; % INFILE4 = [DIR_DATA_ORIG, 'PCAWG/signatures/sigProfiler/PCAWG_sigProfiler_',signatureType,'_signatures_in_samples.csv'];
    INFILE5 = sProperties.GENES_DRIVERS_PCAWG_S3;                                   % [DIR_DATA_ORIG, 'genes/TableS3_panorama_driver_mutations_ICGC_samples.controlled.txt'];
end

saveFileData = ['save/dataPCAWG_',expressionType,'.mat'];
if (runAgain || ~exist(saveFileData, 'file'))
    tic
    %%
    fprintf('Computing %s...\n', saveFileData);
    %% Samples
    tableSampleSheet = readtable(INFILE1);
    tableSampleSheet_tumourWGS = tableSampleSheet(strcmp(tableSampleSheet.library_strategy, 'WGS') & contains(tableSampleSheet.dcc_specimen_type, 'tumour'),:);
    tableSampleSheet_normalWGS = tableSampleSheet(strcmp(tableSampleSheet.library_strategy, 'WGS') & contains(tableSampleSheet.dcc_specimen_type, 'Normal'),:);
    tableSampleSheet_tumourRNA = tableSampleSheet(strcmp(tableSampleSheet.library_strategy, 'RNA-Seq') & contains(tableSampleSheet.dcc_specimen_type, 'tumour'),:);
    tableSampleSheet_normalRNA = tableSampleSheet(strcmp(tableSampleSheet.library_strategy, 'RNA-Seq') & ~contains(tableSampleSheet.dcc_specimen_type, 'tumour'),:);
    tableSamplesPCAWG = unique(tableSampleSheet_tumourWGS(:, {'icgc_sample_id', 'icgc_donor_id', 'icgc_specimen_id', 'dcc_project_code', 'dcc_specimen_type', 'aliquot_id', 'donor_wgs_exclusion_white_gray'}));
    tableSamplesPCAWG.aliquot_id_WGS_tumour = tableSamplesPCAWG.aliquot_id; tableSamplesPCAWG.aliquot_id = [];
    nSamplesPCAWG = size(tableSamplesPCAWG, 1);
    tableSamplesPCAWG.iSampleAll = (1:nSamplesPCAWG)';
    %% Mapping normal WGS aliquot_id_WGS_normal
    [isOK2, indexIntoSheet] = ismember(tableSamplesPCAWG.icgc_donor_id, tableSampleSheet_normalWGS.icgc_donor_id);
    fprintf('%d samples total, %d samples in normal WGS, %d samples not | using icgc_donor_id.\n', nSamplesPCAWG, sum(isOK2), sum(~isOK2)); % 2955 samples total, 2955 samples in normal WGS, 0 samples not | using icgc_donor_id.
    tableSamplesPCAWG.aliquot_id_WGS_normal = cell(nSamplesPCAWG, 1); tableSamplesPCAWG.aliquot_id_WGS_normal(:) = {''};
    tableSamplesPCAWG.aliquot_id_WGS_normal(isOK2) = tableSampleSheet_normalWGS.aliquot_id(indexIntoSheet(isOK2));    
    %% Mapping RNA-seq aliquot_id_RNA: higher priority is to have the same icgc_specimen_id, but if not available (in 227 WGS samples), we use a different RNA sample with the same icgc_donor_id
    %     [~, tableSampleSheet_tumourRNA.iSampleAll] = ismember(tableSampleSheet_tumourRNA.icgc_donor_id, tableSamples.icgc_donor_id);
    %     tableSamples.nRNA = histcounts(tableSampleSheet_tumourRNA.iSampleAll, 1:nSamples+1)'; % Always the first sample of icgc_donor_id is used
    %     fprintf('%d donors have more than 1 RNA sample\n', sum(tableSamples.nRNA>1)); % 17 donors have more than 1 RNA sample
    %
    [isOK2, indexIntoSheet] = ismember(tableSamplesPCAWG.icgc_donor_id, tableSampleSheet_tumourRNA.icgc_donor_id);
    fprintf('%d samples total, %d samples in RNA-seq, %d samples not | using icgc_donor_id.\n', nSamplesPCAWG, sum(isOK2), sum(~isOK2)); % 2955 samples total, 1286 samples in RNA-seq, 1669 samples not | using icgc_donor_id.
    tableSamplesPCAWG.has_RNA = isOK2;
    tableSamplesPCAWG.aliquot_id_RNA = cell(nSamplesPCAWG, 1); tableSamplesPCAWG.aliquot_id_RNA(:) = {''};
    tableSamplesPCAWG.aliquot_id_RNA(isOK2) = tableSampleSheet_tumourRNA.aliquot_id(indexIntoSheet(isOK2));
    tableSamplesPCAWG.has_RNA_specimenMatched = false(size(tableSamplesPCAWG, 1), 1);
    
    [isOK2, indexIntoSheet] = ismember(tableSamplesPCAWG.icgc_specimen_id, tableSampleSheet_tumourRNA.icgc_specimen_id);
    fprintf('%d samples total, %d samples in RNA-seq, %d samples not | using icgc_specimen_id.\n', nSamplesPCAWG, sum(isOK2), sum(~isOK2)); % 2955 samples total, 1059 samples in RNA-seq, 1896 samples not | using icgc_specimen_id.
    tableSamplesPCAWG.aliquot_id_RNA(isOK2) = tableSampleSheet_tumourRNA.aliquot_id(indexIntoSheet(isOK2));
    tableSamplesPCAWG.has_RNA_specimenMatched(isOK2) = true;
    %% Expression
    tableExpression = readtable(INFILE2);
    tableExpressionSamples = table();
    tableExpressionSamples.aliquot_id_RNA = cellfun(@(x) strrep(x(end-35:end), '_', '-'), tableExpression.Properties.VariableNames(2:end)', 'UniformOutput', false); % Those starting with 0 had x added at the beginning --> which we need to remove
    % Just checking which samples these are...
    isOK = ismember(tableExpressionSamples.aliquot_id_RNA, tableSampleSheet_tumourRNA.aliquot_id);
    fprintf('%d samples with expression found in tableSampleSheet_tumourRNA, %d NOT.\n', sum(isOK), sum(~isOK)); % 1305 samples with expression found in tableSampleSheet_tumourRNA, 216 NOT.
    isOK = ismember(tableExpressionSamples.aliquot_id_RNA, tableSampleSheet_normalRNA.aliquot_id);
    fprintf('%d samples with expression found in tableSampleSheet_normalRNA, %d NOT.\n', sum(isOK), sum(~isOK)); % 161 samples with expression found in tableSampleSheet_normalRNA, 1360 NOT.
    isOK = ismember(tableExpressionSamples.aliquot_id_RNA, tableSampleSheet_tumourRNA.aliquot_id) | ismember(tableExpressionSamples.aliquot_id_RNA, tableSampleSheet_normalRNA.aliquot_id);
    fprintf('%d samples with expression found in tableSampleSheet_tumourRNA or tableSampleSheet_normalRNA, %d NOT.\n', sum(isOK), sum(~isOK)); % 1466 samples with expression found in tableSampleSheet_tumourRNA or tableSampleSheet_normalRNA, 55 NOT. - WHAT ARE THESE 55 SAMPLES?
    %%
    [tableExpressionSamples.isInOurSamples, tableExpressionSamples.iSampleAll] = ismember(tableExpressionSamples.aliquot_id_RNA, tableSamplesPCAWG.aliquot_id_RNA);
    fprintf('%d samples with expression found in our sample list (of WGS tumours), %d NOT.\n', sum(tableExpressionSamples.isInOurSamples), sum(~tableExpressionSamples.isInOurSamples)); % 1284 samples with expression found in our sample list, 237 NOT.
    matExpression = table2array(tableExpression(:, [false;tableExpressionSamples.isInOurSamples]'));
    lstExpressionGenes = tableExpression.feature;
    tableExpressionSamples = tableExpressionSamples(tableExpressionSamples.isInOurSamples,:);
    matExpressionGenesSamplesPCAWG = NaN*ones(length(lstExpressionGenes), nSamplesPCAWG);
    matExpressionGenesSamplesPCAWG(:, tableExpressionSamples.iSampleAll) = matExpression;
    tableSamplesPCAWG.isInExpressionTable = false(nSamplesPCAWG, 1); tableSamplesPCAWG.isInExpressionTable(tableExpressionSamples.iSampleAll) = true;
    clear matExpression tableExpression tableExpressionSamples
    clear isOK2 indexIntoSheet tableSampleSheet tableSampleSheet_tumourWGS tableSampleSheet_tumourRNA
    %% Gene names
    tableGencodeGenes = readtable(INFILE3);
    tableGencodeGenes.Properties.VariableNames = {'chr', 'pos0', 'pos1', 'geneSymbol', 'geneType1', 'strand', 'geneNameGencode', 'geneType'};
    tableGencodeGenes.geneType1 = [];
    tableGencodeGenes.geneNameGencodeFirstPart = cellfun(@(x) x(1:min(strfind(x, '.'))-1), tableGencodeGenes.geneNameGencode, 'UniformOutput', false);
    geneNameGencode = lstExpressionGenes;
    tableExpressionGenes = table(geneNameGencode);
    [isOK, index] = ismember(tableExpressionGenes.geneNameGencode, tableGencodeGenes.geneNameGencode);
    fprintf('%d expression genes found in gencode, %d NOT.\n', sum(isOK), sum(~isOK)); % 57820 expression genes found in gencode, 0 NOT.
    if (sum(~isOK)>0), error('ERROR: some genes not found.\n'); end
    tableExpressionGenes = tableGencodeGenes(index, :); nGenesPCAWG = size(tableExpressionGenes, 1);
    clear tableGencodeGenes lstExpressionGenes geneNameGencode isOK index
    %% Driver mutations    
    tableDriverMutationsPCAWG = readtable(INFILE5);
    [~, tableDriverMutationsPCAWG.iSampleAllPCAWG] = ismember(tableDriverMutationsPCAWG.sample_id, tableSamplesPCAWG.aliquot_id_WGS_tumour);
    for iRow = 1:size(tableDriverMutationsPCAWG, 1)
        iSPRow = tableDriverMutationsPCAWG.iSampleAllPCAWG(iRow);
        tableDriverMutationsPCAWG.icgc_sample_id{iRow} = tableSamplesPCAWG.icgc_sample_id{iSPRow};
        tableDriverMutationsPCAWG.icgc_donor_id{iRow} = tableSamplesPCAWG.icgc_donor_id{iSPRow};
        tableDriverMutationsPCAWG.sampleName{iRow} = ['PCAWG_', strrep(tableSamplesPCAWG.dcc_project_code{iSPRow}, '-', '_'), '_', tableSamplesPCAWG.icgc_donor_id{iSPRow}, '_', tableSamplesPCAWG.icgc_sample_id{iSPRow}];
        tableDriverMutationsPCAWG.donorName{iRow} = ['PCAWG_', strrep(tableSamplesPCAWG.dcc_project_code{iSPRow}, '-', '_'), '_', tableSamplesPCAWG.icgc_donor_id{iSPRow}];
    end
    %tableDriverMutationsPCAWG = tableDriverMutationsPCAWG(~strcmp(tableDriverMutationsPCAWG.pos, 'x'), :);
    tableDriverMutationsPCAWG.posNumeric = cellfun(@str2double, tableDriverMutationsPCAWG.pos);
    chrPrenumeric = tableDriverMutationsPCAWG.chr; 
    chrPrenumeric(strcmp(chrPrenumeric, 'X')) = {'23'}; 
    chrPrenumeric(strcmp(chrPrenumeric, 'Y')) = {'24'};
    tableDriverMutationsPCAWG.chrNumeric = cellfun(@str2double, chrPrenumeric);
    %% Mutational signatures
    lstSignatureTypes = {'SBS', 'ID', 'DBS'};
    for iSignatureType = 1:length(lstSignatureTypes)
        signatureType = lstSignatureTypes{iSignatureType};
        INFILE4 = [INDIR4, 'PCAWG_sigProfiler_',signatureType,'_signatures_in_samples.csv']; % INFILE4 = [DIR_DATA_ORIG, 'PCAWG/signatures/sigProfiler/PCAWG_sigProfiler_',signatureType,'_signatures_in_samples.csv'];
        tableMS = readtable(INFILE4);
        nSignatures = size(tableMS, 2)-3;
        [tableMS.isInOurSamples, tableMS.iSampleAll] = ismember(tableMS.SampleNames, tableSamplesPCAWG.icgc_specimen_id);
        tableMS = tableMS(tableMS.isInOurSamples,:);
        lstCols = 3+(1:nSignatures);
        sSignatures.(signatureType).matSamplesSignatures = NaN*ones(nSamplesPCAWG, nSignatures);
        sSignatures.(signatureType).matSamplesSignatures(tableMS.iSampleAll, :) = table2array(tableMS(:,lstCols));
        sSignatures.(signatureType).signatureNames = tableMS.Properties.VariableNames(lstCols);
        sSignatures.(signatureType).nSignatures = nSignatures;
        tableSamplesPCAWG.(['hasMS_',(signatureType)]) = false(nSamplesPCAWG, 1); tableSamplesPCAWG.hasMS_SBS(tableMS.iSampleAll) = true;
        tableSamplesPCAWG.(['accuracyMS_',(signatureType)]) = NaN*ones(nSamplesPCAWG, 1); tableSamplesPCAWG.accuracyMS_SBS(tableMS.iSampleAll) = tableMS.Accuracy;
    end
    sSignatures.lstSignatureTypes = [lstSignatureTypes, {'together'}];
    sSignatures.together.matSamplesSignatures = [sSignatures.SBS.matSamplesSignatures, sSignatures.ID.matSamplesSignatures, sSignatures.DBS.matSamplesSignatures];
    sSignatures.together.signatureNames = [sSignatures.SBS.signatureNames, sSignatures.ID.signatureNames, sSignatures.DBS.signatureNames];
    sSignatures.together.nSignatures = sSignatures.SBS.nSignatures + sSignatures.ID.nSignatures + sSignatures.DBS.nSignatures;
    clear tableMS iSignatureType signatureType nSignatures lstCols
    %%
    %     grpstats(tableSamples(:,{'donor_wgs_exclusion_white_gray', 'has_RNA'}), 'donor_wgs_exclusion_white_gray', 'sum')
    %     %                  donor_wgs_exclusion_white_gray    GroupCount    sum_has_RNA
    %     %     Whitelist            {'Whitelist'}                2703          1190
    %     %     Graylist             {'Graylist' }                  75            32
    %     %     Excluded             {'Excluded' }                 177            64
    %     tmp = tableSamples(ismember(tableSamples.donor_wgs_exclusion_white_gray,{'Excluded'}),:);
    %% 121 non-unique samples
    tableDonorsPCAWG = grpstats(tableSamplesPCAWG(:, {'icgc_donor_id', 'dcc_project_code', 'has_RNA', 'isInExpressionTable'}), {'icgc_donor_id', 'dcc_project_code'}, 'sum');
    tableDonorsPCAWG.has_RNA = tableDonorsPCAWG.sum_has_RNA>0;
    %tableDonorsPCAWG(tableDonorsPCAWG.GroupCount>1,:)
    nDonors = size(tableDonorsPCAWG, 1);
    tableDonorsPCAWG.iDonor = (1:nDonors)';
    [~, tableSamplesPCAWG.iDonor] = ismember(tableSamplesPCAWG.icgc_donor_id, tableDonorsPCAWG.icgc_donor_id);
    tableDonorsPCAWG.nSamplesTumourWGS = histcounts(tableSamplesPCAWG.iDonor, 1:nDonors+1)'; 
    %% Selecting unique samples priorities (highest to lowest): has_RNA_specimenMatched > has_RNA > Primary in dcc_specimen_type > anything
    tableDonorsPCAWG.iSampleAll = NaN*tableDonorsPCAWG.iDonor;
    tableDonorsPCAWG.iSampleAll(tableSamplesPCAWG.iDonor) = tableSamplesPCAWG.iSampleAll;
    isOK = contains(tableSamplesPCAWG.dcc_specimen_type, 'Primary');
    tableDonorsPCAWG.iSampleAll(tableSamplesPCAWG.iDonor(isOK)) = tableSamplesPCAWG.iSampleAll(isOK);
    isOK = tableSamplesPCAWG.has_RNA;
    tableDonorsPCAWG.iSampleAll(tableSamplesPCAWG.iDonor(isOK)) = tableSamplesPCAWG.iSampleAll(isOK);
    isOK = tableSamplesPCAWG.has_RNA_specimenMatched;
    tableDonorsPCAWG.iSampleAll(tableSamplesPCAWG.iDonor(isOK)) = tableSamplesPCAWG.iSampleAll(isOK);
    tableSamplesPCAWG.isUsed = false(nSamplesPCAWG, 1);
    tableSamplesPCAWG.isUsed(tableDonorsPCAWG.iSampleAll) = true;
    %     tmp = tableSamples(ismember(tableSamples.icgc_donor_id, tableDonors.icgc_donor_id(tableDonors.GroupCount>1)),:);
    %     tmp = sortrows(tmp,'icgc_donor_id','ascend');
    fprintf('%d donors, %d have multiple samples, %d samples in total in these donors\n', nDonors, sum(tableDonorsPCAWG.nSamplesTumourWGS>1), sum(tableDonorsPCAWG.nSamplesTumourWGS(tableDonorsPCAWG.nSamplesTumourWGS>1))); % 2834 donors, 63 have multiple samples, 184 samples in total in these donors
    %%
    toc
    save(saveFileData, 'tableSamplesPCAWG', 'tableDonorsPCAWG', 'nSamplesPCAWG', 'tableExpressionGenes', 'matExpressionGenesSamplesPCAWG', 'sSignatures', 'nGenesPCAWG', 'tableDriverMutationsPCAWG');
else
    fprintf('Loading data from %s...\n', saveFileData);
    load(saveFileData, 'tableSamplesPCAWG', 'tableDonorsPCAWG', 'nSamplesPCAWG', 'tableExpressionGenes', 'matExpressionGenesSamplesPCAWG', 'sSignatures', 'nGenesPCAWG', 'tableDriverMutationsPCAWG');
end
%%
% 
% tableProjects = grpstats(tableDonorsPCAWG(:, {'dcc_project_code', 'has_RNA', 'GroupCount'}), {'dcc_project_code'}, 'sum');
% tableProjects = sortrows(tableProjects,'dcc_project_code','ascend');
% tmp = tableProjects(tableProjects.sum_has_RNA>0,:);
