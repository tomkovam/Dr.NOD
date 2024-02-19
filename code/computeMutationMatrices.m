function [matGenesSamplesNMut_SNVs_highCADD, matGenesSamplesNMut_SNVs, matGenesSamplesNMut_INDEL, matUESamplesIsMut_SNVs_highCADD, matUESamplesIsMut_INDEL, matUESamplesIsMut_SNVs] = ...
    computeMutationMatrices(runAgain, suffix, minCADD_PHRED, tableMutations, matMutationsEnhancers, nSamples, tableUniqueEnhancers, matUniqueEnhancersGenes, doSave, sProperties)
%% Annotates enhancers with counts of samples that have a mutation in the enhancer
% Note that matGenesSamplesNMut_SNVs is the number of SNV-mutated UE per every gene and sample 
% (not the total number of all mutations in the regulatory space - this will differ in cases when one regulatory space has multiple mutations in the same sample!)

fileNameMutationMatrices = [sProperties.DIRECTORY_SAVE, '/mutationsMatrices/mutationsMatrices_', suffix, '_', num2str(minCADD_PHRED),'.mat'];
if (~runAgain && exist(fileNameMutationMatrices, 'file'))
    fprintf('Loading a %s...\n', fileNameMutationMatrices);
    load(fileNameMutationMatrices, 'matGenesSamplesNMut_SNVs_highCADD', 'matGenesSamplesNMut_SNVs', 'matGenesSamplesNMut_INDEL', 'matUESamplesIsMut_SNVs_highCADD', 'matUESamplesIsMut_INDEL', 'matUESamplesIsMut_SNVs');
else
    t1 = tic;
    fprintf('annotateEnhancersByMutations...\n');
    %%
    nUE = size(tableUniqueEnhancers, 1);
    % nSamples = size(tableSamples, 1);
    % nGenes = size(matUniqueEnhancersGenes, 2);
    %%
    matUESamplesIsMut_SNVs_highCADD = false(nUE, nSamples);
    matUESamplesIsMut_INDEL = false(nUE, nSamples);
    matUESamplesIsMut_SNVs = false(nUE, nSamples);
    
    for iUE = 1:nUE                                                                                                                                                         % We take the mutations of that iUE unique enhancer and look at their sample list
        matUESamplesIsMut_SNVs_highCADD(iUE, tableMutations.iSample(~tableMutations.isIndel & tableMutations.isHighCADD & matMutationsEnhancers(:, iUE)==1)) = true;        % high-CADD SNVs
        matUESamplesIsMut_SNVs(iUE, tableMutations.iSample(~tableMutations.isIndel & matMutationsEnhancers(:, iUE)==1)) = true;                                             % all-CADD SNVs
        matUESamplesIsMut_INDEL(iUE, tableMutations.iSample(tableMutations.isIndel & matMutationsEnhancers(:, iUE)==1)) = true;                                             % indels
    end
    toc(t1)
    %%
    %     tableUniqueEnhancers.nMutSamples_SNVs = sum(matUESamplesIsMut_SNVs(:,~tableSamples.isExcluded), 2);
    %     tableUniqueEnhancers.nMutSamples_SNVs_highCADD = sum(matUESamplesIsMut_SNVs_highCADD(:,~tableSamples.isExcluded), 2);
    %     tableUniqueEnhancers.nMutSamples_INDELs = sum(matUESamplesIsMut_INDEL(:,~tableSamples.isExcluded), 2);
    %% Genes x samples: how many UE of the given gene are mutated in the given sample? This is faster then the older approach (and I have checked it gives identical results)
    matGenesSamplesNMut_SNVs_highCADD = matUniqueEnhancersGenes'*matUESamplesIsMut_SNVs_highCADD; % matGenesSamplesNMut_SNVs_highCADD(iGene,iSample) ~ how many UE regulating iGene are mutated in sample iSample --> then we care if at least one
    matGenesSamplesNMut_SNVs          = matUniqueEnhancersGenes'*matUESamplesIsMut_SNVs;
    matGenesSamplesNMut_INDEL         = matUniqueEnhancersGenes'*matUESamplesIsMut_INDEL; 
    %%
    toc(t1)
    if (doSave)
        fprintf('Saving %s...\n', fileNameMutationMatrices);
        createDir(fileparts(fileNameMutationMatrices));
        save(fileNameMutationMatrices, 'matGenesSamplesNMut_SNVs_highCADD', 'matGenesSamplesNMut_SNVs', 'matGenesSamplesNMut_INDEL', 'matUESamplesIsMut_SNVs_highCADD', 'matUESamplesIsMut_INDEL', 'matUESamplesIsMut_SNVs');
    else
        fprintf('NOT saving %s...\n', fileNameMutationMatrices);
    end
end