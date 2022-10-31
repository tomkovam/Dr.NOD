function [tableDataForBMM, tableGenes_annotations] = prepare_tableDataForBMM(nMutSamplesInEnhancersPerGene, tableGenes_annotations, tableGenes_mean_trinucleotides)

tableGenes_annotations.mutationFrequency = nMutSamplesInEnhancersPerGene./tableGenes_annotations.nTheoreticalMutations_PHRED_geqCUTOFF;
tableDataForBMM = [tableGenes_mean_trinucleotides, tableGenes_annotations(:,{'mean_mfInFlanks', 'nPositionsInEnhancers', 'mean_GC', 'mean_replicationTiming', 'mean_baseActivity', 'mutationFrequency'})]; % The last column is the response variable

if (ismember("mean_replicationTiming", tableDataForBMM.Properties.VariableNames)) % The zero means a missing value in this data set
    tableDataForBMM.mean_replicationTiming(tableDataForBMM.mean_replicationTiming==0) = NaN;
end

tableDataForBMM.mutationFrequency(isinf(tableDataForBMM.mutationFrequency)) = NaN; % This can happen when nTheoreticalMutations_PHRED_geqCUTOFF=0 and nMutSamplesInEnhancersPerGene>0 because of indels (should not happen for highCADD SNVs)
