function [tableGenes_pValues, stats] = computePValuePerGene(runAgain, suffix, minCADD_PHRED, exclusionType, matExpressionGenesSamples, matGenesSamplesNMut_SNVs_highCADD, matGenesSamplesNMut_INDEL, matCNV_genesSamples, ...
    matUniqueEnhancersGenes, tableGenesNasserExpressed, tableGenes_annotations, tableGenes_mean_trinucleotdies, tableSamples, tableUniqueEnhancers, verbose, sProperties)
% Computes p-value that enhancers of that gene are being more mutated than expected and p-value of expression being different between samples with a
% mutation in one of the enhancers vs WT.

fileNamePValuePerGene = [sProperties.DIRECTORY_SAVE, '/pValuePerGene/pValuePerGene_', suffix, '_', num2str(minCADD_PHRED), '_', exclusionType, '.mat']; 
if (~runAgain && exist(fileNamePValuePerGene, 'file'))
    fprintf('Loading %s...\n', fileNamePValuePerGene);
    load(fileNamePValuePerGene, 'tableGenes_pValues', 'stats');
else
    if (~exist('matExpressionGenesSamples', 'var'))
        error('ERROR: matExpressionGenesSamples parameter is missing');
    end
    t1 = tic;
    fprintf('computePValuePerGene %s...\n', suffix);
    
    if (~isequal(tableGenesNasserExpressed.geneName, tableGenes_annotations.geneName)), error('ERROR gene names in tableGenesNasserExpressed and tableGenes_annotations do not match.'); end
    if (~isequal(size(tableGenesNasserExpressed, 1), size(tableGenes_mean_trinucleotdies, 1))), error('ERROR number of rows of tableGenesNasserExpressed and tableGenes_mean_trinucleotdies do not match.'); end
    
    s1 = warning('error', 'stats:LinearModel:RankDefDesignMat'); s2 = warning('error', 'stats:glmfit:IterationLimit'); s3 = warning('error', 'stats:glmfit:BadScaling');
    
    lstMutTypes = {'SNVs_highCADD', 'SNVs_highCADD_INDELs'};  nTypes = length(lstMutTypes); % OLD: 'SNVs', 'INDELs', 'SNVs_INDELs',
    lstSelTypesMuts = [1,2];

    tableGenes_pValues = tableGenesNasserExpressed(:,'geneName');
    nGenes = size(tableGenes_pValues, 1);
    nUsedSamples = sum(~tableSamples.isExcluded);
    lstExcluded = find(tableSamples.isExcluded)';
    
    matGenesSamplesIsMut_SNVs_highCADD = matGenesSamplesNMut_SNVs_highCADD>0; % Genes x samples: how many UE of the given gene are mutated in the given sample? This is faster then the older approach (and I have checked it gives identical results)
    matGenesSamplesIsMut_INDEL = matGenesSamplesNMut_INDEL>0;

    tmp.SNVs_highCADD = matGenesSamplesIsMut_SNVs_highCADD;
    tmp.SNVs_highCADD_INDELs = matGenesSamplesIsMut_SNVs_highCADD | matGenesSamplesIsMut_INDEL;

    tableGenes_pValues.nMutSamplesInEnhancers_SNVs_highCADD = sum(tmp.SNVs_highCADD(:,~tableSamples.isExcluded), 2); % I checked that this gives the same results as the previous for-loop computation
    tableGenes_pValues.nMutSamplesInEnhancers_SNVs_highCADD_INDELs = sum(tmp.SNVs_highCADD_INDELs(:,~tableSamples.isExcluded), 2);
    
    tableGenes_pValues.nMutSamplesInEnhancers_hasRNA_SNVs_highCADD = sum(tmp.SNVs_highCADD(:,~tableSamples.isExcluded & tableSamples.has_RNA), 2);
    tableGenes_pValues.nMutSamplesInEnhancers_hasRNA_SNVs_highCADD_INDELs = sum(tmp.SNVs_highCADD_INDELs(:,~tableSamples.isExcluded & tableSamples.has_RNA), 2);

    %     sum(sum(matGenesSamplesNMut_SNVs_highCADD(:,~tableSamples.isExcluded)>0))
    %     sum(tableGenes_pValues.nMutSamplesInEnhancers_SNVs_highCADD)        % 1886 correct vs 17412 wrong --> 17335
    
    tableGenes_pValues.nPositionsInEnhancers = matUniqueEnhancersGenes'*tableUniqueEnhancers.nPositions; % Number of positions in all unique enhancers regulating the given gene (I have checked it gives the same results as the previous slower implementation)
    
    for iType = 1:nTypes
        typeName = lstMutTypes{iType};
        lstTests = {'M_basic', 'M_CADDnormalized', 'M_fullModel', 'E_ranksum', 'E_normalLog2', 'E_poisson'};
        for iTest = 1:length(lstTests)
            testName = lstTests{iTest};
            tableGenes_pValues.(['p',testName,'_', typeName]) = NaN*ones(nGenes, 1); % p-value
            tableGenes_pValues.(['e',testName,'_', typeName]) = NaN*ones(nGenes, 1); % estimate (such as estimate in the GLM model or medianDifferenceMutMinusWT for ranksum test)
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
            %             if (iType == 1)
            %                 isMutatedSample = matGenesSamplesIsMut_SNVs_highCADD(iGene,:); %   OLD: isMutatedSample = matGenesSamplesIsMut_SNVs(iGene,:); isMutatedSample = matGenesSamplesIsMut_INDEL(iGene,:);
            %             elseif (iType == 2)
            %                 isMutatedSample = matGenesSamplesIsMut_SNVs_highCADD(iGene,:) | matGenesSamplesIsMut_INDEL(iGene,:); % OLD: isMutatedSample = matGenesSamplesIsMut_SNVs(iGene,:) | matGenesSamplesIsMut_INDEL(iGene,:);
            %             end
            typeName = lstMutTypes{iType};
            isMutatedSample = tmp.(typeName)(iGene,:)' + 0;
            isMutatedSample(lstExcluded) = NaN; % To ignore the excluded samples in this analysis
            %tableGenes_pValues.(['nMutSamplesInEnhancers2_', typeName])(iGene) = sum(isMutatedSample == iGroupMutatedInEnhancer);
            isMutatedSample(isnan(expressionPerSample)) = NaN;
            %tableGenes_pValues.(['nMutSamplesInEnhancers2_hasRNA_', typeName])(iGene) = sum(isMutatedSample == iGroupMutatedInEnhancer);
            if (sum(isMutatedSample==iGroupWT)>0 && sum(isMutatedSample==iGroupMutatedInEnhancer)>0) % There is at least one sample in each group
                tableGenes_pValues.(['pE_ranksum_', typeName])(iGene) = ranksum(expressionPerSample(isMutatedSample==iGroupWT), expressionPerSample(isMutatedSample==iGroupMutatedInEnhancer));
                tableGenes_pValues.(['eE_ranksum_', typeName])(iGene) = median(expressionPerSample(isMutatedSample==iGroupMutatedInEnhancer)) - median(expressionPerSample(isMutatedSample==iGroupWT));

                if (sum(isMutatedSample==iGroupMutatedInEnhancer & ~isnan(expressionPerSample) & ~isnan(CNVperSample))>1)
                    try
                        mdl =  fitglm([isMutatedSample, CNVperSample], log2(1+expressionPerSample));
                        tableGenes_pValues.(['pE_normalLog2_', typeName])(iGene) = mdl.Coefficients.pValue(2); % pValue order: (1) intercept, (2) sampleGroup, (3) CNVperSample
                        tableGenes_pValues.(['eE_normalLog2_', typeName])(iGene) = mdl.Coefficients.Estimate(2);
                    catch
                        %warning('WARNING: problem in computeStatsPerGene GLM_normal_log2expression_  %s', tableGenesResults.geneName{iGene});
                    end
                    try % s = warning('error', 'stats:LinearModel:RankDefDesignMat'); s = warning('error', 'stats:glmfit:IterationLimit'); s = warning('error', 'stats:glmfit:BadScaling'); 'stats:glmfit:BadScaling' [warnMsg, warnId] = lastwarn()
                        mdl =  fitglm([isMutatedSample, CNVperSample], expressionPerSample, 'linear', 'Distribution', 'poisson', 'DispersionFlag', true); % Poisson with overdispersed count variable, as in https://www.nature.com/articles/s41467-019-13929-1 (quasi-Poisson family GLM)
                        tableGenes_pValues.(['pE_poisson_', typeName])(iGene) = mdl.Coefficients.pValue(2); % pValue order: (1) intercept, (2) sampleGroup, (3) CNVperSample
                        tableGenes_pValues.(['eE_poisson_', typeName])(iGene) = mdl.Coefficients.Estimate(2); % This was newly added and needs to be re-run
                    catch
                        %warning('WARNING: problem in computeStatsPerGene GLM_poisson_expression %s', tableGenesResults.geneName{iGene});
                    end
                end
            end
        end
    end
    %%
    stats = struct();
    for iType = lstSelTypesMuts
        typeName = lstMutTypes{iType};
        nMutSamplesInEnhancersPerGene = tableGenes_pValues.(['nMutSamplesInEnhancers_', typeName]);

        % Simpler versions of the mutagenesis model, which we don't use in the end, and therefore do not compute to speed-it-up
        %         averageMutationFrequency = sum(nMutSamplesInEnhancersPerGene)/sum(tableGenes_pValues.nPositionsInEnhancers);
        %         tableGenes_pValues.(['pM_basic_', typeName]) = myBinomTestRightSided(nMutSamplesInEnhancersPerGene, tableGenes_pValues.nPositionsInEnhancers, averageMutationFrequency, 'one');
        %         tableGenes_pValues.(['eM_basic_', typeName]) = log2(nMutSamplesInEnhancersPerGene./(tableGenes_annotations.nTheoreticalMutations_PHRED_geqCUTOFF*averageMutationFrequency));
        %
        %         averageMutationFrequency = sum(nMutSamplesInEnhancersPerGene)/sum(tableGenes_annotations.nTheoreticalMutations_PHRED_geqCUTOFF); % Here we normalize by the number of theoretical mutations with CADD PHRED >= CUTOFF (in all positions in the enhancers regulating given gene)
        %         tableGenes_pValues.(['pM_CADDnormalized_', typeName]) = myBinomTestRightSided(nMutSamplesInEnhancersPerGene, tableGenes_annotations.nTheoreticalMutations_PHRED_geqCUTOFF, averageMutationFrequency, 'one');
        %         tableGenes_pValues.(['eM_CADDnormalized_', typeName]) = log2(nMutSamplesInEnhancersPerGene./(tableGenes_annotations.nTheoreticalMutations_PHRED_geqCUTOFF*averageMutationFrequency));
        
        [xPValues, statsOne] = compute_pM_GLM_geneLevel(nMutSamplesInEnhancersPerGene, tableGenes_annotations, tableGenes_mean_trinucleotdies, nUsedSamples);
        tableGenes_pValues.(['pM_fullModel_', typeName]) = xPValues;
        tableGenes_pValues.(['eM_fullModel_', typeName]) = statsOne.foldChange;
        stats.(typeName) = statsOne;
    end
    %%
    toc(t1)
    createDir(fileparts(fileNamePValuePerGene));
    save(fileNamePValuePerGene, 'tableGenes_pValues', 'stats');
end
%%
if (exist('tableGenesNasserExpressed', 'var'))
    if (~isequal(tableGenesNasserExpressed.geneName, tableGenes_pValues.geneName)), error('ERROR geneNames do not match.'); end
end
