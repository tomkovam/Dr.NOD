function [tableMutations, tableSamples, tableGenesNasserExpressed, tableGencodeGenes, tableCGC, tableDriverMutationsPCAWG, tableDriverGenesPCAWG, ...
    tableEnhancers, tableUniqueEnhancers, tableUniqueEnhancers_regions, matUniqueEnhancersGenes, ...
    matExpressionGenesSamples, matSamplesSignatures, signatureNames, matMutationsEnhancers, matCNV_genesSamples, matSV_genesSamples_minDistance, matSV_genesSamples_nameSVs] = ...
    loadCancerData(runAgain, tissueName, biosampleABC, enhancerAnalysis, doSave, verbose, tissueNameSV, sProperties)

if (~ismember(enhancerAnalysis, {'All', 'Noncoding', 'Slop250bpAll', 'Slop250bpNoncoding'}))
    error('Allowed values for enhancerAnalysis are: All or Noncoding, while %s was used.\n', enhancerAnalysis);
end

expressionType = 'fpkm_uq'; % This is fixed atm, as I did not want to have the save names too long
suffix = [tissueName, '_', biosampleABC, '_', enhancerAnalysis];
saveFileData = ['save/dataCancer_',suffix,'.mat'];
if (runAgain || ~exist(saveFileData, 'file'))
    t1 = tic;
    %%
    fprintf('\n=================\nComputing %s...\n', saveFileData);
    [tableSamplesPCAWG, ~, ~, tableExpressionGenes, matExpressionGenesSamplesPCAWG, sSignatures, tableDriverMutationsPCAWG] = loadDataPCAWG(runAgain, expressionType, sProperties);
    [tableMutations, tableSamples] = loadMutations(runAgain, sProperties, tissueName, tableSamplesPCAWG, tableDriverMutationsPCAWG);
    [tableGenesNasser, tableGencodeGenes, tableChrSizes, tableCGC, tableDriverGenesPCAWG] = loadGenes(runAgain, sProperties, tableDriverMutationsPCAWG);
    [tableEnhancers, tableGenesNasserExpressed, tableUniqueEnhancers, tableUniqueEnhancers_regions, matUniqueEnhancersGenes] = loadDataEnhancers(biosampleABC, enhancerAnalysis, tableGenesNasser, verbose, sProperties);
    %% Gene expression
    [tableGenesNasserExpressed.isInExpression, tableGenesNasserExpressed.indexExpression] = ismember(tableGenesNasserExpressed.geneNameGencode, tableExpressionGenes.geneNameGencode);
    if (min(tableGenesNasserExpressed.isInExpression)<1), error('ERROR: some genes not found.'); end
    if (~isequal(tableSamplesPCAWG.icgc_sample_id(tableSamples.iSamplePCAWG), tableSamples.icgc_sample_id)), error('ERROR: samples do not match.'); end
    matExpressionGenesSamples = matExpressionGenesSamplesPCAWG(tableGenesNasserExpressed.indexExpression, tableSamples.iSamplePCAWG);
    %% Mutational signatures
    signatureNames = sSignatures.together.signatureNames;
    matSamplesSignatures = sSignatures.together.matSamplesSignatures(tableSamples.iSamplePCAWG,:);
    isSignaturePOLE = ismember(signatureNames, {'SBS10a', 'SBS10b', 'DBS3'});
    isSignatureMSI = ismember(signatureNames, {'SBS6', 'SBS14', 'SBS15', 'SBS20', 'SBS21', 'SBS44', 'DBS7', 'DBS10', 'ID7'});
    tableSamples.totalSignatureExposure = sum(matSamplesSignatures, 2);
    tableSamples.percSignaturePOLE_MSI = 100*sum(matSamplesSignatures(:,isSignaturePOLE|isSignatureMSI), 2)./tableSamples.totalSignatureExposure;
    clear matExpressionGenesSamplesPCAWG sSignatures nGenesPCAWG
    %% Intersect mutations and enhancers
    [tableMutations, matMutationsEnhancers] = mutationsInEnhancers(tableMutations, tableUniqueEnhancers_regions, tableChrSizes);
    %%
    doSave_CNV_SV = false; % To save space
    [matCNV_genesSamples, sample_has_CNV] = loadPCAWG_CNV_oneProject(runAgain, suffix, tissueNameSV, biosampleABC, tableGenesNasserExpressed, tableGencodeGenes, tableSamples, sProperties, doSave_CNV_SV);
    [matSV_genesSamples_minDistance, matSV_genesSamples_nameSVs, sample_has_SV] = loadPCAWG_SV_oneProject(runAgain, suffix, tissueNameSV, biosampleABC, tableGenesNasserExpressed, tableGencodeGenes, tableSamples, sProperties, doSave_CNV_SV);
    tableSamples.has_CNV = sample_has_CNV;
    tableSamples.has_SV = sample_has_SV;
    tableSamples.project = cellfun(@(x) x(7:strfind(x,'_DO')-1), tableSamples.sampleName, 'UniformOutput', false);
    %% Are any of our enhancer mutations annotated as PCAWG driver mutations?
    if (verbose)
        isOK1 = tableMutations.isPCAWGDriver & tableMutations.iUniqueEnhancer>0;
        fprintf('%s %s: PCAWG DRIVER MUTATIONS in tableMutations: %d mutations, %d SNVs in %s enhancers.\n', tissueName, biosampleABC, sum(tableMutations.isPCAWGDriver), sum(isOK1), enhancerAnalysis);
        if (sum(isOK1)>0)
            fprintf('Known driver mutations in enhancer regions and tableSNVs:\n');
            tmp1 = tableMutations(isOK1,:)
            tmp2 = tableDriverMutationsPCAWG(ismember(tableDriverMutationsPCAWG.posNumeric, tmp1.pos1) & ismember(tableDriverMutationsPCAWG.chrNumeric, tmp1.chrNumeric),:)
        end
    end
    %% Save
    toc(t1)
    if (doSave)
        save(saveFileData, 'tableMutations', 'tableSamples', 'tableGenesNasser', 'tableGencodeGenes', 'tableCGC', 'tableDriverMutationsPCAWG', 'tableDriverGenesPCAWG', ...
        'tableEnhancers', 'tableGenesNasserExpressed', 'tableUniqueEnhancers', 'tableUniqueEnhancers_regions', 'matUniqueEnhancersGenes', ...
        'matExpressionGenesSamples', 'matSamplesSignatures', 'signatureNames', 'matMutationsEnhancers', 'matCNV_genesSamples', 'matSV_genesSamples_minDistance', 'matSV_genesSamples_nameSVs');
    end
else
    fprintf('Loading data from %s...\n', saveFileData);
    load(saveFileData, 'tableMutations', 'tableSamples', 'tableGenesNasser', 'tableGencodeGenes', 'tableCGC', 'tableDriverMutationsPCAWG', 'tableDriverGenesPCAWG', ...
        'tableEnhancers', 'tableGenesNasserExpressed', 'tableUniqueEnhancers', 'tableUniqueEnhancers_regions', 'matUniqueEnhancersGenes', ...
        'matExpressionGenesSamples', 'matSamplesSignatures', 'signatureNames', 'matMutationsEnhancers', 'matCNV_genesSamples', 'matSV_genesSamples_minDistance', 'matSV_genesSamples_nameSVs');
end