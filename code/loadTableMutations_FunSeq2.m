function [tableMutations_FunSeq2] = loadTableMutations_FunSeq2(runAgain, tissueName, biosampleABC, tableMutations, tableSamples, sProperties)
%% Annotates the SNV mutations with FunSeq2 predictions of altered TFBS

enhancerAnalysis = sProperties.enhancerAnalysis;
minCADD_PHRED = sProperties.minCADD_PHRED;
exclusionType = sProperties.exclusionType;

saveFileData = ['save/tableMutations_FunSeq2/tableMutations_FunSeq2_',tissueName,'_', biosampleABC, '_', enhancerAnalysis, '_', num2str(minCADD_PHRED), '_', exclusionType, '.mat'];
if (runAgain || ~exist(saveFileData, 'file'))
    t0 = tic;
    %% 
    fprintf('Computing %s...\n', saveFileData);

    lstProjects = unique(tableSamples.projectName); 


    tableMutations.chr = strrep(cellstr(num2str(tableMutations.chrNumeric,'chr%d')), ' ', '');
    tableMutations.chr(tableMutations.chrNumeric == 23) = {'chrX'};
    tableMutations.chr(tableMutations.chrNumeric == 24) = {'chrY'};
    tableMutations.mutID = strcat(tableMutations.chr, '_', strrep(cellstr(num2str(tableMutations.pos1)), ' ', ''), '_', tableMutations.ref, '_', tableMutations.alt);

    tableMutationsAnnotatedFunSeq2 = tableMutations;
    tableMutationsAnnotatedFunSeq2.FunSeq2_isAnnotated = false(size(tableMutationsAnnotatedFunSeq2, 1), 1);
    tableMutationsAnnotatedFunSeq2.FunSeq2_gerp = NaN*ones(size(tableMutationsAnnotatedFunSeq2, 1), 1);
    tableMutationsAnnotatedFunSeq2.FunSeq2_noncodingScore = NaN*ones(size(tableMutationsAnnotatedFunSeq2, 1), 1);
    tableMutationsAnnotatedFunSeq2.FunSeq2_isMOTIFBR = false(size(tableMutationsAnnotatedFunSeq2, 1), 1);
    tableMutationsAnnotatedFunSeq2.FunSeq2_isMOTIFG = false(size(tableMutationsAnnotatedFunSeq2, 1), 1);
    
    % Example row:
    % chr, pos0, pos1, ref, alt, gerp, cds (Yes/No), motif, noncodingScore
    % chr1	861250	861251	G	A	-3.16	No	MOTIFBR=CTCF,DHS,SUZ12#BHLHE40_disc2#861236#861253#-#3#0.050633#0.451477,MOTIFG=AHR::ARNT_3#861249#861255#-#5#7.208#2.813	1.15963101894413   

    lstCols1 = {'chr', 'pos0', 'pos1', 'ref', 'alt', 'FunSeq2_gerp', 'FunSeq2_noncodingScore', 'FunSeq2_isMOTIFBR', 'FunSeq2_isMOTIFG'};

    for iProject = 1:length(lstProjects)
        projectName = lstProjects{iProject};
        t1 = tic;
        try
            tableAnnotatedFunSeq2 = readtable([sProperties.FUNSEQ2_DIR, 'intersectedB.',projectName,'.motif.FunSeq2.bed.txt'], 'ReadVariableNames', true); % ['data/FunSeq2/annotatedPCAWG/intersectedB.',projectName,'.motif.FunSeq2.bed.txt']
            tableAnnotatedFunSeq2.Properties.VariableNames = {'chr', 'pos0', 'pos1', 'ref', 'alt', 'FunSeq2_chr', 'FunSeq2_pos0', 'FunSeq2_pos1', 'FunSeq2_ref', 'FunSeq2_alt', 'FunSeq2_gerp', 'FunSeq2_cds', 'FunSeq2_motifAnalysis', 'FunSeq2_noncodingScore'};
            %tableAnnotatedFunSeq2.FunSeq2_isAnnotated = 1 + 0*tableAnnotatedFunSeq2.pos0;
            tableAnnotatedFunSeq2.FunSeq2_cds = []; %strcmp(tmp.FunSeq2_cds, 'Yes');
            tableAnnotatedFunSeq2.FunSeq2_isMOTIFBR = contains(tableAnnotatedFunSeq2.FunSeq2_motifAnalysis, 'MOTIFBR');
            tableAnnotatedFunSeq2.FunSeq2_isMOTIFG = contains(tableAnnotatedFunSeq2.FunSeq2_motifAnalysis, 'MOTIFG');
            tableAnnotatedFunSeq2.FunSeq2_motifAnalysis = [];
            isOK = arrayfun(@(iRow) strcmp(tableAnnotatedFunSeq2.ref{iRow}, tableAnnotatedFunSeq2.FunSeq2_ref{iRow}) & contains(tableAnnotatedFunSeq2.FunSeq2_alt{iRow}, tableAnnotatedFunSeq2.alt{iRow}), (1:size(tableAnnotatedFunSeq2, 1))', 'UniformOutput', true);
            tableAnnotatedFunSeq2 = tableAnnotatedFunSeq2(isOK,lstCols1);
        catch
            warning('FunSeq2 file for project %s not found.', projectName);
            continue
        end

        if (~isempty(tableAnnotatedFunSeq2))
            tableAnnotatedFunSeq2 = unique(tableAnnotatedFunSeq2);

            tableAnnotatedFunSeq2_unique = grpstats(unique(tableAnnotatedFunSeq2), {'chr', 'pos0', 'pos1', 'ref', 'alt'}, 'max');
            fprintf('%d rows with multiple different FunSeq2 annotations.\n', sum(tableAnnotatedFunSeq2_unique.GroupCount>1));
            %tableAnnotatedFunSeq2_unique(tableAnnotatedFunSeq2_unique.GroupCount>1,:)
            %tableAnnotatedFunSeq2(ismember(tableAnnotatedFunSeq2.pos1, tableAnnotatedFunSeq2_unique.pos1(tableAnnotatedFunSeq2_unique.GroupCount>1)),:)
            tableAnnotatedFunSeq2_unique.GroupCount = [];
            tableAnnotatedFunSeq2_unique.Properties.VariableNames = lstCols1; % {'chr', 'pos0', 'pos1', 'ref', 'alt', 'FunSeq2_gerp', 'FunSeq2_noncodingScore', 'FunSeq2_isMOTIFBR', 'FunSeq2_isMOTIFG'};
            tableAnnotatedFunSeq2_unique.chrNumeric = cellfun(@(x) str2double(x(4:end)), tableAnnotatedFunSeq2_unique.chr);
            tableAnnotatedFunSeq2_unique.chrNumeric(strcmp(tableAnnotatedFunSeq2_unique.chr, 'chrX')) = 23;
            tableAnnotatedFunSeq2_unique.chrNumeric(strcmp(tableAnnotatedFunSeq2_unique.chr, 'chrY')) = 24;
            tableAnnotatedFunSeq2_unique.mutID = strcat(tableAnnotatedFunSeq2_unique.chr, '_', strrep(cellstr(num2str(tableAnnotatedFunSeq2_unique.pos1)), ' ', ''), '_', tableAnnotatedFunSeq2_unique.ref, '_', tableAnnotatedFunSeq2_unique.alt);

            [isInFunSeq2, indexFunSeq2] = ismember(tableMutationsAnnotatedFunSeq2.mutID, tableAnnotatedFunSeq2_unique.mutID);
            lstCols2 = {'FunSeq2_gerp', 'FunSeq2_noncodingScore', 'FunSeq2_isMOTIFBR', 'FunSeq2_isMOTIFG'};
            tableMutationsAnnotatedFunSeq2.FunSeq2_isAnnotated(isInFunSeq2) = 1;
            tableMutationsAnnotatedFunSeq2(isInFunSeq2,lstCols2) = tableAnnotatedFunSeq2_unique(indexFunSeq2(isInFunSeq2), lstCols2);

            isMutation = ~tableMutationsAnnotatedFunSeq2.isIndel & ~tableMutationsAnnotatedFunSeq2.isExcluded & ismember(tableMutationsAnnotatedFunSeq2.iSample, find(strcmp(tableSamples.projectName, projectName)));
            fprintf('%s: %.1f%% mutations have FunSeq2 annotation.\n', projectName, 100*mean(tableMutationsAnnotatedFunSeq2.FunSeq2_isAnnotated(isMutation)));
            toc(t1)
        end
    end

    fprintf('%s: %.1f%% mutations have FunSeq2 annotation.\n', tissueName, 100*mean(tableMutationsAnnotatedFunSeq2.FunSeq2_isAnnotated));

    if (size(tableMutations, 1) ~= size(tableMutationsAnnotatedFunSeq2, 1)) % This shouldn't happen
        error('Rows do not match.');
    end

    tableAllEnhancers = readtable(sProperties.ABC_ENHANCER_GENE_MAPS_MERGED); % 'data/Nasser2021/Slop250bpAllPredictions.together.merged.bed.txt'
    tableAllEnhancers.Properties.VariableNames = {'chr', 'pos0', 'pos1'};
    tableAllEnhancers.chrNumeric = cellfun(@(x) str2double(x(4:end)), tableAllEnhancers.chr);
    tableAllEnhancers.chrNumeric(strcmp(tableAllEnhancers.chr, 'chrX')) = 23;
    tableAllEnhancers.chrNumeric(strcmp(tableAllEnhancers.chr, 'chrY')) = 24;

    tableGenesNasser = readtable(sProperties.GENES_ABC, 'delimiter', '\t'); % 'data/genes/Nasser2021_genesTSS.formatted.txt'
    tableGenesNasser.Properties.VariableNames = {'nEnhancersAllBiosamples', 'chr', 'geneName', 'TSS'};
    tableGenesNasser = tableGenesNasser(:,{'chr', 'geneName', 'TSS', 'nEnhancersAllBiosamples'});
    tableGenesNasser.chrNumeric = cellfun(@(x) str2double(x(4:end)), tableGenesNasser.chr);
    tableGenesNasser.chrNumeric(strcmp(tableGenesNasser.chr, 'chrX')) = 23;
    tableGenesNasser.chrNumeric(strcmp(tableGenesNasser.chr, 'chrY')) = 24;

    nMutations = size(tableMutationsAnnotatedFunSeq2, 1);
    tableMutationsAnnotatedFunSeq2.isInABCEnhancer = true(nMutations, 1);
    tableMutationsAnnotatedFunSeq2.isNearTSS_250bp = false(nMutations, 1);
    tableMutationsAnnotatedFunSeq2.isInABCEnhancer_thisTissue = false(nMutations, 1);
    tableMutationsAnnotatedFunSeq2.isNearTSS_250bp_thisTissue = false(nMutations, 1);

    for iMutation = 1:nMutations
        isOK = tableAllEnhancers.chrNumeric == tableMutationsAnnotatedFunSeq2.chrNumeric(iMutation) & tableAllEnhancers.pos0 <= tableMutationsAnnotatedFunSeq2.pos1(iMutation) & tableAllEnhancers.pos1 >= tableMutationsAnnotatedFunSeq2.pos0(iMutation);
        tableMutationsAnnotatedFunSeq2.isInABCEnhancer(iMutation) = sum(isOK)>0;
        isOK = tableGenesNasser.chrNumeric == tableMutationsAnnotatedFunSeq2.chrNumeric(iMutation) & abs(tableGenesNasser.TSS - tableMutationsAnnotatedFunSeq2.pos1(iMutation)) <= 250;
        tableMutationsAnnotatedFunSeq2.isNearTSS_250bp(iMutation) = sum(isOK)>0;
    end
    %%
    lstCols3 = {'pos1', 'FunSeq2_isAnnotated', 'FunSeq2_isMOTIFBR', 'FunSeq2_isMOTIFG', 'FunSeq2_gerp', 'FunSeq2_noncodingScore', 'isInABCEnhancer', 'isNearTSS_250bp'};
    tableMutations_FunSeq2 = tableMutationsAnnotatedFunSeq2(:,lstCols3);
    %%
    toc(t0)
    createDir(fileparts(saveFileData));
    save(saveFileData, 'tableMutations_FunSeq2');
else
    fprintf('Loading data from %s...\n', saveFileData);
    load(saveFileData, 'tableMutations_FunSeq2');
end

if (~isequal(tableMutations.pos1, tableMutations_FunSeq2.pos1)), error('Rows in tableMutations and tableMutations_FunSeq2 do not match.'); end
