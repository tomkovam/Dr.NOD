function [pM, stats] = compute_pM_GLM_geneLevel(nMutSamplesInEnhancersPerGene, tableGenes_annotations, tableGenes_mean_trinucleotides, nUsedSamples, univariableCutoff, verbose)
%% Computes the pM p-values (for scoreM) for each gene, i.e., whether the regulatory space of the gene is more mutated than expected (taking high-CADD mutations into consideration only)

if (~exist('univariableCutoff', 'var'))
    univariableCutoff = 0.001;
end
if (~exist('verbose', 'var'))
    verbose = false;
end

tableGenes_annotations.mutationFrequency = nMutSamplesInEnhancersPerGene./tableGenes_annotations.nTheoreticalMutations_PHRED_geqCUTOFF;
tableDataForBMM = [tableGenes_mean_trinucleotides, tableGenes_annotations(:,{'mean_mfInFlanks', 'nPositionsInEnhancers', 'mean_GC', 'mean_replicationTiming', 'mean_baseActivity', 'mutationFrequency'})]; % The last column is the response variable


if (ismember("mean_replicationTiming", tableDataForBMM.Properties.VariableNames)) % The zero means a missing value in this data set
    tableDataForBMM.mean_replicationTiming(tableDataForBMM.mean_replicationTiming==0) = NaN;
end

tableDataForBMM.mutationFrequency(isinf(tableDataForBMM.mutationFrequency)) = NaN; % This can happen when nTheoreticalMutations_PHRED_geqCUTOFF=0 and nMutSamplesInEnhancersPerGene>0 because of indels (should not happen for highCADD SNVs)

tablePredictors = table();
tablePredictors.predictor = tableDataForBMM.Properties.VariableNames';
nPredictors = size(tablePredictors, 1);

tablePredictors.univariablePValue = NaN*ones(nPredictors, 1);
for iCol = 1:nPredictors
    try
        mdl =  fitglm(tableDataForBMM(:,[iCol,end]), 'linear', 'Distribution', 'poisson', 'DispersionFlag', true);
        tablePredictors.univariablePValue(iCol) = mdl.coefTest;
    catch
        warning('PROBLEM in column %s', tablePredictors.predictor{iCol});
    end
end

try
    mdl =  fitglm(tableDataForBMM(:,tablePredictors.univariablePValue<univariableCutoff), 'linear', 'Distribution', 'poisson', 'DispersionFlag', true);
    if (verbose)
        tablePredictors(tablePredictors.univariablePValue<univariableCutoff,:)
        mdl
        fprintf('Used predictors: %d\nExplained deviance: %.g\n', sum(tablePredictors.univariablePValue<univariableCutoff), mdl.Rsquared.Deviance);
    end
    expected_mf = predict(mdl, tableDataForBMM(:,1:end-1));

    %% BinomTest(n, k, 1 - (1-p)^s, ’one’)
    n = nMutSamplesInEnhancersPerGene;
    k = nUsedSamples;
    s = tableGenes_annotations.nTheoreticalMutations_PHRED_geqCUTOFF;
    p = expected_mf/nUsedSamples;
    
    pM = myBinomTestRightSided(n, k, 1 - ((1-p).^s), 'one'); 
    pM(pM>1) = 1; % In the case of NaN inputs
    %%
    warning('off', 'MATLAB:nearlySingularMatrix') % So that mdl.coefTest runs without warnings
    stats.explainedDevaince = mdl.Rsquared.Deviance;
    stats.explainedVariance = mdl.Rsquared.Ordinary;
    stats.coefTest = mdl.coefTest;
    stats.tablePredictors = tablePredictors;
    stats.expected_mf = expected_mf;
    stats.observed_mf = tableGenes_annotations.mutationFrequency;
    stats.foldChange = (stats.observed_mf./stats.expected_mf);
    stats.log2FC = log2(stats.observed_mf./stats.expected_mf);
    stats.univariableCutoff = univariableCutoff;
    %stats.mdl = mdl; % Too large
catch
    warning('Badly scaled data or some other issues.');
end
