function [tableGenes_pValues_hyperUE, stats_hyperUE, tableUE_annotations_hyperUE, tmpUE, statsUE, matUESamplesIsMut_SNVs_highCADD_hyperUE, matUESamplesIsMut_SNVs_highCADD_INDEL_hyperUE] = computePValuePerGene_hypermutatedEnhancersOnly(runAgain, suffix, minCADD_PHRED, exclusionType, matExpressionGenesSamples, matCNV_genesSamples, ...
    matUniqueEnhancersGenes, tableGenesNasserExpressed, tableGenes_annotations, tableGenes_mean_trinucleotdies, matUESamplesIsMut_SNVs_highCADD, matUESamplesIsMut_INDEL, tableUE_annotations, tableUE_mean_trinucleotdies, tableSamples, verbose)
%% Recomputes the y-axis (expression) p-values with only samples mutated in hypermutated enhancers (by a fullModel of background mutagenesis, computer per enhancer)
% Used only for plotting the genome-view examples

alpha_hyperUE = 0.05;

fileNamePValuePerGene = ['save/pValuePerGene/pValuePerGene_hyperUE_', suffix, '_', num2str(minCADD_PHRED), '_', exclusionType, '_', num2str(alpha_hyperUE), '.mat'];
if (~runAgain && exist(fileNamePValuePerGene, 'file'))
    fprintf('Loading %s...\n', fileNamePValuePerGene);
    load(fileNamePValuePerGene, 'tableGenes_pValues_hyperUE', 'stats_hyperUE', 'tableUE_annotations_hyperUE', 'tmpUE', 'statsUE', 'matUESamplesIsMut_SNVs_highCADD_hyperUE', 'matUESamplesIsMut_SNVs_highCADD_INDEL_hyperUE');
else
    if (~exist('matExpressionGenesSamples', 'var'))
        error('ERROR: matExpressionGenesSamples parameter is missing');
    end
    t1 = tic;
    fprintf('computePValuePerGene %s...\n', suffix);
    
    if (~isequal(tableGenesNasserExpressed.geneName, tableGenes_annotations.geneName)), error('ERROR gene names in tableGenesNasserExpressed and tableGenes_annotations do not match.'); end
    if (~isequal(size(tableGenesNasserExpressed, 1), size(tableGenes_mean_trinucleotdies, 1))), error('ERROR number of rows of tableGenesNasserExpressed and tableGenes_mean_trinucleotdies do not match.'); end
    
    s1 = warning('error', 'stats:LinearModel:RankDefDesignMat'); s2 = warning('error', 'stats:glmfit:IterationLimit'); s3 = warning('error', 'stats:glmfit:BadScaling');
    
    tableGenes_pValues_hyperUE = tableGenesNasserExpressed(:,'geneName');
    nGenes = size(tableGenes_pValues_hyperUE, 1);
    nUsedSamples = sum(~tableSamples.isExcluded);
    lstExcluded = find(tableSamples.isExcluded)';
    
    lstMutTypes = {'SNVs_highCADD', 'SNVs_highCADD_INDELs'};  nTypes = length(lstMutTypes); % OLD: 'SNVs', 'INDELs', 'SNVs_INDELs',
    lstSelTypesMuts = [1,2];

    clear tmpUE statsUE
    tmpUE = struct();
    statsUE = struct();
    tmpUE.SNVs_highCADD = matUESamplesIsMut_SNVs_highCADD;
    tmpUE.SNVs_highCADD_INDELs = matUESamplesIsMut_SNVs_highCADD | matUESamplesIsMut_INDEL;

    tableUE_annotations_hyperUE = tableUE_annotations;

    tableUE_annotations_hyperUE.mean_mfInFlanks = tableUE_annotations_hyperUE.mfInFlanks;
    tableUE_annotations_hyperUE.nPositionsInEnhancers = tableUE_annotations_hyperUE.nPositions;


    tableUE_annotations_hyperUE.nMutSamplesInEnhancers_SNVs_highCADD = sum(tmpUE.SNVs_highCADD(:,~tableSamples.isExcluded), 2); % I checked that this gives the same results as the previous for-loop computation
    tableUE_annotations_hyperUE.nMutSamplesInEnhancers_SNVs_highCADD_INDELs = sum(tmpUE.SNVs_highCADD_INDELs(:,~tableSamples.isExcluded), 2);

    tableUE_annotations_hyperUE.nMutSamplesInEnhancers_hasRNA_SNVs_highCADD = sum(tmpUE.SNVs_highCADD(:,~tableSamples.isExcluded & tableSamples.has_RNA), 2);
    tableUE_annotations_hyperUE.nMutSamplesInEnhancers_hasRNA_SNVs_highCADD_INDELs = sum(tmpUE.SNVs_highCADD_INDELs(:,~tableSamples.isExcluded & tableSamples.has_RNA), 2);

    for iType = lstSelTypesMuts
        typeName = lstMutTypes{iType};
        nMutSamplesPerUE = tableUE_annotations_hyperUE.(['nMutSamplesInEnhancers_', typeName]);

        averageMutationFrequency = sum(nMutSamplesPerUE)/sum(tableUE_annotations_hyperUE.nPositions);
        tableUE_annotations_hyperUE.(['pM_basic_', typeName]) = myBinomTestRightSided(nMutSamplesPerUE, tableUE_annotations_hyperUE.nPositions, averageMutationFrequency, 'one');
        tableUE_annotations_hyperUE.(['eM_basic_', typeName]) = log2(nMutSamplesPerUE./(tableUE_annotations_hyperUE.nTheoreticalMutations_PHRED_geqCUTOFF*averageMutationFrequency));

        averageMutationFrequency = sum(nMutSamplesPerUE)/sum(tableUE_annotations_hyperUE.nTheoreticalMutations_PHRED_geqCUTOFF); % Here we normalize by the number of theoretical mutations with CADD PHRED >= CUTOFF (in all positions in the enhancers regulating given gene)
        tableUE_annotations_hyperUE.(['pM_CADDnormalized_', typeName]) = myBinomTestRightSided(nMutSamplesPerUE, tableUE_annotations_hyperUE.nTheoreticalMutations_PHRED_geqCUTOFF, averageMutationFrequency, 'one');
        tableUE_annotations_hyperUE.(['eM_CADDnormalized_', typeName]) = log2(nMutSamplesPerUE./(tableUE_annotations_hyperUE.nTheoreticalMutations_PHRED_geqCUTOFF*averageMutationFrequency));

        [xPValues, statsOne] = compute_pM_GLM_geneLevel(nMutSamplesPerUE, tableUE_annotations_hyperUE, tableUE_mean_trinucleotdies, nUsedSamples);
        tableUE_annotations_hyperUE.(['pM_fullModel_', typeName]) = xPValues;
        tableUE_annotations_hyperUE.(['eM_fullModel_', typeName]) = statsOne.foldChange;
        statsUE.(typeName) = statsOne;
    end

    matUESamplesIsMut_SNVs_highCADD_hyperUE = matUESamplesIsMut_SNVs_highCADD;
    matUESamplesIsMut_SNVs_highCADD_INDEL_hyperUE = matUESamplesIsMut_SNVs_highCADD | matUESamplesIsMut_INDEL;

    matUESamplesIsMut_SNVs_highCADD_hyperUE(~(tableUE_annotations_hyperUE.pM_fullModel_SNVs_highCADD<alpha_hyperUE),:) = false;
    matUESamplesIsMut_SNVs_highCADD_INDEL_hyperUE(~(tableUE_annotations_hyperUE.pM_fullModel_SNVs_highCADD_INDELs<alpha_hyperUE),:) = false;

    clear tmp stats_hyperUE
    tmp.SNVs_highCADD = (matUniqueEnhancersGenes'*matUESamplesIsMut_SNVs_highCADD_hyperUE)>0;
    tmp.SNVs_highCADD_INDELs = (matUniqueEnhancersGenes'*matUESamplesIsMut_SNVs_highCADD_INDEL_hyperUE)>0;
    stats_hyperUE = struct();

    for iType = 1:nTypes
        typeName = lstMutTypes{iType};
        lstTests = {'M_basic', 'M_CADDnormalized', 'M_fullModel', 'E_ranksum', 'E_normalLog2', 'E_poisson'};
        for iTest = 1:length(lstTests)
            testName = lstTests{iTest};
            tableGenes_pValues_hyperUE.(['p',testName,'_', typeName]) = NaN*ones(nGenes, 1); % p-value
            tableGenes_pValues_hyperUE.(['e',testName,'_', typeName]) = NaN*ones(nGenes, 1); % estimate (such as estimate in the GLM model or medianDifferenceMutMinusWT for ranksum test)
        end
    end

    for iGene = 1:nGenes
        if (verbose)
            if (mod(iGene, 1e3)==1), fprintf('Gene %d...\n', iGene); end
        end
        expressionPerSample = matExpressionGenesSamples(iGene, :)';
        expressionPerSample(lstExcluded) = NaN; % To ignore the excluded samples in this analysis
        CNVperSample=matCNV_genesSamples(iGene, :)';

        for iType = lstSelTypesMuts
            iGroupWT = 0;
            iGroupMutatedInEnhancer = 1;
            typeName = lstMutTypes{iType};
            isMutatedSample = tmp.(typeName)(iGene,:)' + 0;
            isMutatedSample(lstExcluded) = NaN; % To ignore the excluded samples in this analysis
            isMutatedSample(isnan(expressionPerSample)) = NaN;
            if (sum(isMutatedSample==iGroupWT)>0 && sum(isMutatedSample==iGroupMutatedInEnhancer)>0) % There is at least one sample in each group
                tableGenes_pValues_hyperUE.(['pE_ranksum_', typeName])(iGene) = ranksum(expressionPerSample(isMutatedSample==iGroupWT), expressionPerSample(isMutatedSample==iGroupMutatedInEnhancer));
                tableGenes_pValues_hyperUE.(['eE_ranksum_', typeName])(iGene) = median(expressionPerSample(isMutatedSample==iGroupMutatedInEnhancer)) - median(expressionPerSample(isMutatedSample==iGroupWT));

                if (sum(isMutatedSample==iGroupMutatedInEnhancer & ~isnan(expressionPerSample) & ~isnan(CNVperSample))>1)
                    try
                        mdl =  fitglm([isMutatedSample, CNVperSample], log2(1+expressionPerSample));
                        tableGenes_pValues_hyperUE.(['pE_normalLog2_', typeName])(iGene) = mdl.Coefficients.pValue(2); % pValue order: (1) intercept, (2) sampleGroup, (3) CNVperSample
                        tableGenes_pValues_hyperUE.(['eE_normalLog2_', typeName])(iGene) = mdl.Coefficients.Estimate(2);
                    catch
                        %warning('WARNING: problem in computeStatsPerGene GLM_normal_log2expression_  %s', tableGenesResults.geneName{iGene});
                    end
                    try % s = warning('error', 'stats:LinearModel:RankDefDesignMat'); s = warning('error', 'stats:glmfit:IterationLimit'); s = warning('error', 'stats:glmfit:BadScaling'); 'stats:glmfit:BadScaling' [warnMsg, warnId] = lastwarn()
                        mdl =  fitglm([isMutatedSample, CNVperSample], expressionPerSample, 'linear', 'Distribution', 'poisson', 'DispersionFlag', true); % Poisson with overdispersed count variable, as in https://www.nature.com/articles/s41467-019-13929-1 (quasi-Poisson family GLM)
                        tableGenes_pValues_hyperUE.(['pE_poisson_', typeName])(iGene) = mdl.Coefficients.pValue(2); % pValue order: (1) intercept, (2) sampleGroup, (3) CNVperSample
                        tableGenes_pValues_hyperUE.(['eE_poisson_', typeName])(iGene) = mdl.Coefficients.Estimate(2); % This was newly added and needs to be re-run
                    catch
                        %warning('WARNING: problem in computeStatsPerGene GLM_poisson_expression %s', tableGenesResults.geneName{iGene});
                    end
                end
            end
        end
    end
    %%
    toc(t1)
    createDir(fileparts(fileNamePValuePerGene));
    save(fileNamePValuePerGene, 'tableGenes_pValues_hyperUE', 'stats_hyperUE', 'tableUE_annotations_hyperUE', 'tmpUE', 'statsUE', 'matUESamplesIsMut_SNVs_highCADD_hyperUE', 'matUESamplesIsMut_SNVs_highCADD_INDEL_hyperUE');
end
%%
if (exist('tableGenesNasserExpressed', 'var'))
    if (~isequal(tableGenesNasserExpressed.geneName, tableGenes_pValues_hyperUE.geneName)), error('ERROR geneNames do not match.'); end
end
