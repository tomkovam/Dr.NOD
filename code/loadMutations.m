function [tableMutations, tableSamples] = loadMutations(runAgain, sProperties, tissueName, tableSamplesPCAWG, tableDriverMutationsPCAWG)

if (~exist('tableSamplesPCAWG', 'var'))
    [tableSamplesPCAWG, ~, ~, ~, ~, ~, tableDriverMutationsPCAWG] = loadDataPCAWG();
end
DIR_INTERSECTED_MUT=sProperties.PCAWG_MUTATIONS_DIR; %[DIR_DATA_NEW, 'out03_inEnhancersPromoters_250bp/'];
%%
saveFileData = ['save/dataMutations_',tissueName,'.mat'];
if (runAgain || ~exist(saveFileData, 'file'))
    tic
    %% Files were copied from cluster, gunzipped using files = gunzip('*.gz');, only the unzipped ones were kept, and .bedlike --> .txt were renamed in TotalCommander.
    fprintf('Computing %s...\n', saveFileData);
    
    tableMutations = table();
    tableSamples = table();
    inFile = sProperties.PCAWG_PROJECT_NAMES; %[DIR_DATA_NEW, 'tableProjects.xlsx']
    tableProjects = readtable(inFile);
    lstProjects = find(strcmp(tableProjects.TISSUE, tissueName))';
    fprintf('%s: %d projects\n', tissueName, length(lstProjects));
    tableProjects(lstProjects,:)

    for iProject = lstProjects
        PROJECT_SAVE_NAME = tableProjects.PROJECT_SAVE_NAME{iProject};
        %% Load tableSamplesOneProject
        inFile = [sProperties.PCAWG_SAMPLES_DIR, 'samples.', PROJECT_SAVE_NAME, '.WGS.txt']; % [DIR_DATA_NEW, 'out04_samples/samples.', PROJECT_SAVE_NAME, '.WGS.txt']
        tableSamplesOneProject = readtable(inFile, 'delimiter', '\t');
        tableSamplesOneProject.Properties.VariableNames = {'iSampleOneProject', 'sampleName'};
        for iSample = 1:size(tableSamplesOneProject)
            tmp = strsplit(tableSamplesOneProject.sampleName{iSample}, '_');
            tableSamplesOneProject.icgc_sample_id{iSample} = tmp{end};
            tableSamplesOneProject.icgc_donor_id{iSample} = tmp{end-1};
        end
        [isOK, tableSamplesOneProject.iSamplePCAWG] = ismember(tableSamplesOneProject.icgc_sample_id, tableSamplesPCAWG.icgc_sample_id);
        if (sum(~isOK)>0), error('ERROR: %d samples in %s not found\n', sum(~isOK), PROJECT_SAVE_NAME); end
        tableSamplesOneProject.isUsedOnePerDonor = tableSamplesPCAWG.isUsed(tableSamplesOneProject.iSamplePCAWG);
        nSamplesOneProject = size(tableSamplesOneProject, 1);
        tableSamplesOneProject.aliquotID_WGS = cell(nSamplesOneProject, 1); tableSamplesOneProject.aliquotID_WGS(:) = {''};
        tableSamplesOneProject.aliquotID_RNA = cell(nSamplesOneProject, 1); tableSamplesOneProject.aliquotID_RNA(:) = {''};
        tableSamplesOneProject.aliquotID_WGS(isOK) = tableSamplesPCAWG.aliquot_id_WGS_tumour(tableSamplesOneProject.iSamplePCAWG(isOK));
        tableSamplesOneProject.aliquotID_RNA(isOK) = tableSamplesPCAWG.aliquot_id_RNA(tableSamplesOneProject.iSamplePCAWG(isOK));
        tableSamplesOneProject.has_RNA = ~strcmp(tableSamplesOneProject.aliquotID_RNA, '');
        fprintf('%s: %d samples. %d samples have PCAWG RNA (%d not).\n', PROJECT_SAVE_NAME, nSamplesOneProject, sum(tableSamplesOneProject.has_RNA), sum(~tableSamplesOneProject.has_RNA));
        tableSamplesOneProject.isTCGA = contains(tableSamplesOneProject.sampleName, '_US_');
        [~, tableDriverMutationsPCAWG.iSample] = ismember(tableDriverMutationsPCAWG.icgc_donor_id, tableSamplesOneProject.icgc_donor_id);
        fprintf('%d rows (%d samples) of tableDriverMutations found in our table (which has %d nonTCGA and %d TCGA samples).\n', ...
            sum(tableDriverMutationsPCAWG.iSample>0), length(unique(tableDriverMutationsPCAWG.iSample(tableDriverMutationsPCAWG.iSample>0))), sum(~tableSamplesOneProject.isTCGA), sum(tableSamplesOneProject.isTCGA));
        %% Load tableSNVs
        tableSNVs = readtable([DIR_INTERSECTED_MUT, 'substitutions.',PROJECT_SAVE_NAME,'.annotated.txt'], 'delimiter', '\t');
        tableSNVs.Properties.VariableNames = {'chr', 'pos0', 'pos1', 'iSample', 'sampleName', 'ref', 'alt', 'supportingReads', 'totalReads', 'CADD_RawScore', 'CADD_PHRED', 'trinucleotide', 'iPattern', 'strand', 'gene', 'context'};
        isOppositeStrand = strcmp(tableSNVs.strand, '-');
        originalREF = tableSNVs.ref;
        originalALT = tableSNVs.alt;
        tableSNVs.ref(isOppositeStrand & strcmp(originalREF, 'C')) = {'G'};
        tableSNVs.ref(isOppositeStrand & strcmp(originalREF, 'T')) = {'A'};
        tableSNVs.alt(isOppositeStrand & strcmp(originalALT, 'A')) = {'T'};
        tableSNVs.alt(isOppositeStrand & strcmp(originalALT, 'C')) = {'G'};
        tableSNVs.alt(isOppositeStrand & strcmp(originalALT, 'G')) = {'C'};
        tableSNVs.alt(isOppositeStrand & strcmp(originalALT, 'T')) = {'A'};
        %fprintf('SNVs converted %s to reverse complement\n', num2sepNumStr(sum(isOppositeStrand)));
        tableSNVs.isIndel = false(size(tableSNVs, 1), 1);
        %% Load tableINDELs and combine them to tableMutationsOneCT
        try
            tableINDELs = readtable([DIR_INTERSECTED_MUT, 'indels.',PROJECT_SAVE_NAME,'.bedlike.txt'], 'delimiter', '\t'); nIndels = size(tableINDELs, 1);
            tableINDELs.Properties.VariableNames = {'chr', 'pos0', 'pos1', 'iSample', 'sampleName', 'ref', 'alt', 'supportingReads', 'totalReads'};

            tableINDELs.CADD_RawScore = NaN*ones(nIndels, 1);
            tableINDELs.CADD_PHRED = NaN*ones(nIndels, 1);
            tableINDELs.trinucleotide = cell(nIndels, 1); tableINDELs.trinucleotide(:) = {''};
            tableINDELs.iPattern = NaN*ones(nIndels, 1);
            tableINDELs.strand = cell(nIndels, 1); tableINDELs.strand(:) = {''};
            tableINDELs.gene = cell(nIndels, 1); tableINDELs.gene(:) = {''};
            tableINDELs.context = cell(nIndels, 1); tableINDELs.context(:) = {''};
            tableINDELs.isIndel = true(nIndels, 1);

            tableMutationsOneCT = [tableSNVs; tableINDELs];
        catch
            error('ERROR: %s issue in indels.\n', PROJECT_SAVE_NAME);
            tableMutationsOneCT = tableSNVs;
        end
        clear tableSNVs tableINDELs
        %% Check that sampleNames and iSample match
        if (~isequal(tableMutationsOneCT.sampleName, tableSamplesOneProject.sampleName(tableMutationsOneCT.iSample))), error('ERROR: sampleNames and iSample do not match.'); end
        %%
        tableMutationsOneCT = unique(tableMutationsOneCT);
        tableMutationsOneCT.isUsedOnePerDonor = ismember(tableMutationsOneCT.iSample, tableSamplesOneProject.iSampleOneProject(tableSamplesOneProject.isUsedOnePerDonor));
        fprintf('%.1f%% mutations belong to isUsedOnePerDonor samples\n', 100*mean(tableMutationsOneCT.isUsedOnePerDonor));
        %% VAF and order by VAF within e
        tableMutationsOneCT.chrNumeric = cellfun(@(x) str2double(x(4:end)), tableMutationsOneCT.chr);
        tableMutationsOneCT.chrNumeric(strcmp(tableMutationsOneCT.chr, 'chrX')) = 23;
        tableMutationsOneCT.chrNumeric(strcmp(tableMutationsOneCT.chr, 'chrY')) = 24;
        tableMutationsOneCT.VAF = tableMutationsOneCT.supportingReads ./ tableMutationsOneCT.totalReads;
        tableMutationsOneCT.iVAF = NaN*tableMutationsOneCT.VAF;
        tableMutationsOneCT.qtlVAF = NaN*tableMutationsOneCT.VAF;
        tableMutationsOneCT = sortrows(tableMutationsOneCT,'iSample','ascend');
        %
        t2 = tic;
        for iSample = 1:max(tableMutationsOneCT.iSample)
            tmpOneSample = tableMutationsOneCT(tableMutationsOneCT.iSample==iSample,:);
            tmpOneSample = sortrows(tmpOneSample,'VAF','descend');
            nMut = size(tmpOneSample,1);
            tmpOneSample.iVAF = (1:nMut)';
            tmpOneSample.qtlVAF = tmpOneSample.iVAF/nMut;
            tableMutationsOneCT(tableMutationsOneCT.iSample==iSample,:) = tmpOneSample;
        end
        toc(t2)
        %%
        tableMutationsOneCT.iSample = NaN*tableMutationsOneCT.iSample;
        tableSamplesOneProject.projectName = cell(size(tableSamplesOneProject, 1), 1); tableSamplesOneProject.projectName(:) = {PROJECT_SAVE_NAME};
        tableMutations = [tableMutations; tableMutationsOneCT];
        tableSamples = [tableSamples; tableSamplesOneProject];
    end
    tableMutations.isIndel = tableMutations.isIndel==1;
    nSamples = size(tableSamples, 1);
    tableSamples.iSample = (1:nSamples)';
    [~, tableMutations.iSample] = ismember(tableMutations.sampleName, tableSamples.sampleName);
    tableSamples.nSNVsIntersected = histcounts(tableMutations.iSample(~tableMutations.isIndel), 1:nSamples+1)';
    tableSamples.nIndelsIntersected = histcounts(tableMutations.iSample(tableMutations.isIndel), 1:nSamples+1)';
    %% Removing redundant columns
    tableMutations.sampleName = [];
    tableMutations.trinucleotide = [];
    tableMutations.CADD_RawScore = [];
    tableMutations.isPlusStrand = strcmp(tableMutations.strand, '+');
    tableMutations.strand = [];
    tableMutations.iVAF = [];
    tableMutations.chr = [];
    %%
    tableMutations.isPCAWGDriver = ismember(tableMutations.chrNumeric, tableDriverMutationsPCAWG.chrNumeric) & ismember(tableMutations.pos1, tableDriverMutationsPCAWG.posNumeric);
    tableMutations.isPCAWGDriver = ismember(tableMutations.chrNumeric, tableDriverMutationsPCAWG.chrNumeric) & ismember(tableMutations.pos1+1, tableDriverMutationsPCAWG.posNumeric) & tableMutations.isIndel;
    tableMutations.driverGene = cell(size(tableMutations, 1), 1); tableMutations.driverGene(:) = {''};
    for iMutation = find(tableMutations.isPCAWGDriver)'
        if (tableMutations.isIndel(iMutation))
            isDriverThisMutation = (tableDriverMutationsPCAWG.chrNumeric == tableMutations.chrNumeric(iMutation)) & (abs(tableDriverMutationsPCAWG.posNumeric-tableMutations.pos1(iMutation))<=1);
        else
            isDriverThisMutation = (tableDriverMutationsPCAWG.chrNumeric == tableMutations.chrNumeric(iMutation)) & (tableDriverMutationsPCAWG.posNumeric == tableMutations.pos1(iMutation));
        end
        tableMutations.driverGene{iMutation} = strjoin(tableDriverMutationsPCAWG.gene(isDriverThisMutation));
    end
    %%
    tableMutations.context = upper(tableMutations.context);
    %%
    toc
    %myPrintMemory
    save(saveFileData, 'tableMutations', 'tableSamples', '-v7.3'); %
else
    fprintf('Loading data from %s...\n', saveFileData);
    load(saveFileData, 'tableMutations', 'tableSamples');
end
fprintf('MUTATIONS %s: %d samples, intersected with enhancers+-250bp and promoters have in median: %s SNVs, %s indels.\n', tissueName, size(tableSamples, 1), ...
    num2sepNumStr(round(median(tableSamples.nSNVsIntersected))), num2sepNumStr(round(median(tableSamples.nIndelsIntersected))));
