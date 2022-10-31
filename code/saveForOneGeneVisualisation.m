function saveForOneGeneVisualisation(tissueName, biosampleABC, geneName, gene_pM, gene_pE, gene_qCombined, tableSamples, matCNV_genesSamples, matExpressionGenesSamples, matGenesSamplesNMut_SNVs_highCADD, ...
    tableMutations, matMutationsEnhancers, iGene, tableGencodeGenes, tableGenesNasserExpressed, matUniqueEnhancersGenes, tableUniqueEnhancers, tableUE_annotations_hyperUE, tableTrinucleotides, exclusionType)

CNVperSample=matCNV_genesSamples(iGene, :)';
expressionPerSample = matExpressionGenesSamples(iGene, :)';
isMutatedSample = matGenesSamplesNMut_SNVs_highCADD(iGene, :)'>0;

if (isempty(expressionPerSample))
    fprintf('No samples with expression for gene %s.\n', geneName);
    return
end

nSamples = size(tableSamples, 1);

sampleGroup = isMutatedSample + 1;
sampleGroup(tableSamples.isExcluded) = NaN;
expressionPerSample(tableSamples.isExcluded) = NaN;

sampleGroupInclWoExpression = sampleGroup;
sampleGroup(isnan(expressionPerSample)) = NaN;

if (sum(sampleGroup==2 & ~isnan(expressionPerSample))==0)
    fprintf('No mutated samples with expression for gene %s (%d mutated samples without expression).\n', geneName, sum(sampleGroup==2));
    return
end
%%
gene_pos0 = tableGencodeGenes.pos0(tableGenesNasserExpressed.iGencode(iGene));
gene_pos1 = tableGencodeGenes.pos1(tableGenesNasserExpressed.iGencode(iGene));
gene_TSS = tableGenesNasserExpressed.TSS(iGene);
gene_strand = tableGenesNasserExpressed.strand_GENCODE(iGene);

lstUniqueEnhancers = find(matUniqueEnhancersGenes(:,iGene)); %unique(tableEnhancers.iUniqueEnhancer(tableEnhancers.iGene == iGeneNasser));
gene_nUEs = length(lstUniqueEnhancers);

tableUniqueEnhancers_oneGene = tableUniqueEnhancers(lstUniqueEnhancers,:);
tableUE_annotations_hyperUE_oneGene = tableUE_annotations_hyperUE(lstUniqueEnhancers,:);


tableUE_annotations_hyperUE_oneGene.foldChangeScoreM = 2.^(tableUE_annotations_hyperUE_oneGene.eM_fullModel_SNVs_highCADD);
tableUE_annotations_hyperUE_oneGene.foldChangeScoreM(isinf(tableUE_annotations_hyperUE_oneGene.eM_fullModel_SNVs_highCADD)) = 0;

isMutPerEnhancer = full(matMutationsEnhancers(:,lstUniqueEnhancers)==1);
isMutOfThisGene = sum(isMutPerEnhancer, 2)>0;

isMutOK = isMutOfThisGene & tableMutations.isHighCADD & ~tableMutations.isExcluded;

tableMutationsThisGene = tableMutations(isMutOK,:);
isMutPerEnhancer = isMutPerEnhancer(isMutOK,:);

if (min(tableMutationsThisGene.iUniqueEnhancer)==0), error('Should be positive'); end
nPatterns = size(tableTrinucleotides, 1);
tableMutationsThisGene.expression = expressionPerSample(tableMutationsThisGene.iSample);
tableMutationsThisGene.yValues = log2(1+tableMutationsThisGene.expression);
tableMutationsThisGene.iPattern(tableMutationsThisGene.isIndel) = nPatterns;
tableMutationsThisGene.patternName = tableTrinucleotides.patternName(tableMutationsThisGene.iPattern);
%%
if (strcmp(exclusionType, 'excludePOLE_MSI'))
    suffix = '';
else
    suffix = ['_', exclusionType];
end
saveFileData = ['save/oneGene/oneGene_', tissueName, '_', biosampleABC, '_', geneName, suffix, '.mat'];
createDir(fileparts(saveFileData));
save(saveFileData, 'gene_pM', 'gene_pE', 'gene_qCombined', 'isMutatedSample', 'expressionPerSample', 'CNVperSample', 'nSamples', 'sampleGroup', 'sampleGroupInclWoExpression', ...
    'gene_pos0', 'gene_pos1', 'gene_TSS', 'gene_strand', 'gene_nUEs', 'tableMutationsThisGene', 'isMutPerEnhancer', 'tableUniqueEnhancers_oneGene', 'tableUE_annotations_hyperUE_oneGene');
