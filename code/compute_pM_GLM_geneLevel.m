function [pM, stats, mdl] = compute_pM_GLM_geneLevel(nMutSamplesInEnhancersPerGene, tableGenes_annotations, tableGenes_mean_trinucleotides, nUsedSamples, univariableCutoff, verbose)
%% Computes the pM p-values (for scoreM) for each gene, i.e., whether the regulatory space of the gene is more mutated than expected (taking high-CADD mutations into consideration only)

if (~exist('univariableCutoff', 'var'))
    univariableCutoff = 0.001;
end
if (~exist('verbose', 'var'))
    verbose = false;
end

[tableDataForBMM, tableGenes_annotations] = prepare_tableDataForBMM(nMutSamplesInEnhancersPerGene, tableGenes_annotations, tableGenes_mean_trinucleotides);

tablePredictors = computeSignificantUnivariablePredictors(tableDataForBMM, univariableCutoff);

try
    mdl =  fitglm(tableDataForBMM(:,tablePredictors.isUsed), 'linear', 'Distribution', 'poisson', 'DispersionFlag', true);
    if (verbose)
        tablePredictors(tablePredictors.isUsed,:)
        mdl
        fprintf('Used predictors: %d\nExplained deviance: %.g\n', sum(tablePredictors.isUsed), mdl.Rsquared.Deviance);
    end
    expected_mf = predict(mdl, tableDataForBMM(:,tablePredictors.isUsed & ~tablePredictors.isResponseVariable));

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
