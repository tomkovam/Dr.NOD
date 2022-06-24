function [pM, stats] = compute_pM_GLM_geneLevel(nMutSamplesInEnhancersPerGene, tableGenes_annotations, tableGenes_mean_trinucleotides, nUsedSamples, univariableCutoff, verbose)

if (~exist('univariableCutoff', 'var'))
    univariableCutoff = 0.001;
end
if (~exist('verbose', 'var'))
    verbose = false;
end
% runEnhancerLevel = false;
% if (runEnhancerLevel)
%     tableUE_annotations.mutationFrequency = tableUniqueEnhancers.nMutSamples_SNVs_highCADD./tableUE_annotations.nTheoreticalMutations_PHRED_geqCUTOFF;
%     tableDataForBMM = [tableUE_mean_trinucleotdies, tableUE_annotations(:,{'mfInFlanks', 'nPositions', 'mean_GC', 'mean_replicationTiming', 'activity_base', 'mutationFrequency'})]; % The last column is the response variable
% else
tableGenes_annotations.mutationFrequency = nMutSamplesInEnhancersPerGene./tableGenes_annotations.nTheoreticalMutations_PHRED_geqCUTOFF;
tableDataForBMM = [tableGenes_mean_trinucleotides, tableGenes_annotations(:,{'mean_mfInFlanks', 'nPositionsInEnhancers', 'mean_GC', 'mean_replicationTiming', 'mean_baseActivity', 'mutationFrequency'})]; % The last column is the response variable
% end

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
% tablePredictors.univariablePValue(end) = 0; % We always need the last column, as that will be the response variable


try
    mdl =  fitglm(tableDataForBMM(:,tablePredictors.univariablePValue<univariableCutoff), 'linear', 'Distribution', 'poisson', 'DispersionFlag', true);
    if (verbose)
        tablePredictors(tablePredictors.univariablePValue<univariableCutoff,:)
        mdl
        fprintf('Used predictors: %d\nExplained deviance: %.g\n', sum(tablePredictors.univariablePValue<univariableCutoff), mdl.Rsquared.Deviance);
    end
    expected_mf = predict(mdl, tableDataForBMM(:,1:end-1));

    % if (runEnhancerLevel)
    %     ypredEnhancerLevel_nMutSamples = ypred.*tableUE_annotations.nTheoreticalMutations_PHRED_geqCUTOFF;
    %     ypred = NaN*tableGenes_annotations.nPositionsInEnhancers;
    %     tmp = NaN*tableGenes_annotations.nPositionsInEnhancers;
    %     tmp2 = NaN*tableGenes_annotations.nPositionsInEnhancers;
    %     for iGene = 1:nGenes
    %         ypred(iGene,:) = sum(ypredEnhancerLevel_nMutSamples(matUniqueEnhancersGenes(:,iGene), :), 1);
    %         tmp(iGene,:) = sum(tableUE_annotations.nTheoreticalMutations_PHRED_geqCUTOFF(matUniqueEnhancersGenes(:,iGene), :), 1);
    %         tmp2(iGene,:) = sum(tableUE_annotations.nPositions(matUniqueEnhancersGenes(:,iGene), :), 1);
    %     end
    %     ypred = ypred./tableGenes_annotations.nTheoreticalMutations_PHRED_geqCUTOFF;
    %     if (max(abs(tableGenes_annotations.nTheoreticalMutations_PHRED_geqCUTOFF-tmp))>0 || max(abs(tableGenes_annotations.nPositionsInEnhancers-tmp2))), error('ERROR: nPositionsInEnhancers or nTheoreticalMutations_PHRED_geqCUTOFF do not match.'); end
    % end
    %%
    % OLD pM = myBinomTestRightSided(nMutSamplesInEnhancersPerGene, tableGenes_annotations.nTheoreticalMutations_PHRED_geqCUTOFF*nUsedSamples, expected_mf/nUsedSamples, 'one'); % myBinomTestRightSided(n, k*s, p, 'one')
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
    %stats.mdl = mdl;
catch
    warning('Badly scaled data or some other issues.');
end
