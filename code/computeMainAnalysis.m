function [tableGenesNasserExpressed, tableGenes_pValues, stats, tableSamples, tableGencodeGenes, ... % levelOutputArguments = 1
    tableGenes_pValues_hyperUE, stats_hyperUE, tableMutations, matMutationsEnhancers, matUniqueEnhancersGenes, matExpressionGenesSamples, ... % levelOutputArguments = 2
    matGenesSamplesNMut_SNVs_highCADD, matGenesSamplesNMut_INDEL, matCNV_genesSamples, tableGenes_annotations, tableGenes_mean_trinucleotdies, tableUniqueEnhancers, ...
    tableDriverMutations, matUESamplesIsMut_SNVs_highCADD, matUESamplesIsMut_INDEL, tableUE_annotations, tableUE_mean_trinucleotdies, matUESamplesIsMut_SNVs_highCADD_hyperUE, tableUE_annotations_hyperUE, matGenesSamplesNMut_SNVs] = ... % levelOutputArguments = 3
    computeMainAnalysis(runAgain, levelOutputArguments, tissueName, biosampleABC, sProperties, tissueNameSV)

% One mutation in tableMutations can be in >1 UE (when neighbouring enhancers are overlapping). Working with matMutationsEnhancers is recommended. This will not cause any double-counting, since we always count max 1 mutation per sample in every computation.

enhancerAnalysis = sProperties.enhancerAnalysis;
minCADD_PHRED = sProperties.minCADD_PHRED;
exclusionType = sProperties.exclusionType;
doSave = sProperties.doSave;
suffix = [tissueName, '_', biosampleABC, '_', enhancerAnalysis];
levelOutputArgumentsOrig = levelOutputArguments;

if (levelOutputArguments <= 2) % To speed it up, we first try to read directly from the precomputed files, assuming that they already exist.
    try
        [tableGenes_pValues, stats] = computePValuePerGene(runAgain, suffix, minCADD_PHRED, exclusionType);
        saveFileData = ['save/dataCancer_',suffix,'.mat'];
        load(saveFileData, 'tableGenesNasserExpressed', 'tableGencodeGenes');
        [tableSamples, tableMutations] = excludeSamples(runAgain, suffix, minCADD_PHRED, exclusionType);
    catch
        levelOutputArguments = 3;
    end
end

if (levelOutputArguments == 2) % To speed it up, we first try to read directly from the precomputed files, assuming that they already exist.
    try
        [tableGenes_pValues_hyperUE, stats_hyperUE] = computePValuePerGene_hypermutatedEnhancersOnly(runAgain, suffix, minCADD_PHRED, exclusionType);
        saveFileData = ['save/dataCancer_',suffix,'.mat'];
        load(saveFileData, 'matMutationsEnhancers', 'matUniqueEnhancersGenes', 'matExpressionGenesSamples');
    catch
        levelOutputArguments = 3;
    end
end


if (levelOutputArguments > 2) % If the files did not exist, or if the user required higher-level of output arguments (meaning more 
    verbose = false;
    fprintf('Loading/computing data for %s %s...\n', tissueName, biosampleABC);
    [tableMutations, tableSamples, tableGenesNasserExpressed, tableGencodeGenes, ~, tableDriverMutations, ~, ...
        ~, tableUniqueEnhancers, tableUniqueEnhancers_regions, matUniqueEnhancersGenes, ...
        matExpressionGenesSamples, ~, ~, matMutationsEnhancers, matCNV_genesSamples] = ...
        loadCancerData(runAgain, tissueName, biosampleABC, enhancerAnalysis, doSave, verbose, tissueNameSV, sProperties);
    %
    [tableSamples, tableMutations] = excludeSamples(runAgain, suffix, minCADD_PHRED, exclusionType, tableSamples, tableMutations, tissueName, biosampleABC, verbose);
    nSamples = size(tableSamples, 1);
    [matGenesSamplesNMut_SNVs_highCADD, matGenesSamplesNMut_SNVs, matGenesSamplesNMut_INDEL, matUESamplesIsMut_SNVs_highCADD, matUESamplesIsMut_INDEL] = ...
        computeMutationMatrices(runAgain, suffix, minCADD_PHRED, tableMutations, matMutationsEnhancers, nSamples, tableUniqueEnhancers, matUniqueEnhancersGenes, doSave);
    %
    [tableGenes_annotations, tableGenes_mean_trinucleotdies, ~, tableUE_annotations, tableUE_mean_trinucleotdies] = ...
        annotateEnhancersByGenomicFeatures(runAgain, suffix, minCADD_PHRED, biosampleABC, enhancerAnalysis, tableSamples, tableUniqueEnhancers, tableGenesNasserExpressed, matUniqueEnhancersGenes, sProperties);
    %
    [tableGenes_pValues, stats] = computePValuePerGene(runAgain, suffix, minCADD_PHRED, exclusionType, matExpressionGenesSamples, matGenesSamplesNMut_SNVs_highCADD, matGenesSamplesNMut_INDEL, matCNV_genesSamples, ...
        matUniqueEnhancersGenes, tableGenesNasserExpressed, tableGenes_annotations, tableGenes_mean_trinucleotdies, tableSamples, tableUniqueEnhancers, verbose);
    %
    tableGenes_pValues_hyperUE = NaN;
    stats_hyperUE = NaN;
    tableUE_annotations_hyperUE = NaN;
    matUESamplesIsMut_SNVs_highCADD_hyperUE = NaN;
    if (levelOutputArgumentsOrig>1)
        [tableGenes_pValues_hyperUE, stats_hyperUE, tableUE_annotations_hyperUE, ~, ~, matUESamplesIsMut_SNVs_highCADD_hyperUE] = ...
            computePValuePerGene_hypermutatedEnhancersOnly(runAgain, suffix, minCADD_PHRED, exclusionType, matExpressionGenesSamples, matCNV_genesSamples, ...
            matUniqueEnhancersGenes, tableGenesNasserExpressed, tableGenes_annotations, tableGenes_mean_trinucleotdies, matUESamplesIsMut_SNVs_highCADD, matUESamplesIsMut_INDEL, tableUE_annotations, tableUE_mean_trinucleotdies, tableSamples, verbose);
    end
    % Add minimal and maximal positions to the unique enhancers:
    tmp = grpstats(tableUniqueEnhancers_regions(:,{'iUniqueEnhancer', 'pos0', 'pos1'}), 'iUniqueEnhancer', {'min', 'max'});
    tableUniqueEnhancers.min_pos0 = NaN*tableUniqueEnhancers.nRegions;
    tableUniqueEnhancers.max_pos1 = NaN*tableUniqueEnhancers.nRegions;
    tableUniqueEnhancers.min_pos0(tmp.iUniqueEnhancer) = tmp.min_pos0;
    tableUniqueEnhancers.max_pos1(tmp.iUniqueEnhancer) = tmp.max_pos1;
    fprintf('\n\n======================== %s %s ========================\n', tissueName, biosampleABC);
end


typeNames = {'SNVs_highCADD', 'SNVs_highCADD_INDELs'};
for iType = 1:length(typeNames)
    typeName = typeNames{iType};
    if (levelOutputArguments <= 2)
        [pM, statsOne] = compute_pValue_pM(runAgain, typeName, tissueName, biosampleABC, enhancerAnalysis, minCADD_PHRED, exclusionType);
    else
        [pM, statsOne] = compute_pValue_pM(runAgain, typeName, tissueName, biosampleABC, enhancerAnalysis, minCADD_PHRED, exclusionType, tableGenes_pValues, tableSamples, tableGenes_annotations, tableGenes_mean_trinucleotdies);
    end
    tableGenes_pValues.(['pM_fullModel_',typeName]) = pM;
    tableGenes_pValues.(['eM_fullModel_', typeName]) = statsOne.foldChange;
    stats.(typeName) = statsOne;
end

