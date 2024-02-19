function [tableSamples, tableMutations] = excludeSamples(runAgain, suffix, minCADD_PHRED, exclusionType, sProperties, tableSamples, tableMutations, tissueName, biosampleABC, verbose)
% Annotated exclued samples (such as POLE-MUT and MSI/MMRd samples) and mutations and annotates high-CADD mutations.

saveFileData = [sProperties.DIRECTORY_SAVE, '/samples/tableSamples_', suffix, '_', num2str(minCADD_PHRED), '_', exclusionType, '.mat'];
if (runAgain || ~exist(saveFileData, 'file'))
    tic
    %%
    fprintf('Computing %s...\n', saveFileData);
    nSamples = size(tableSamples, 1);
    if (strcmp(exclusionType, 'excludePOLE_MSI'))
        tableSamples.isExcluded = tableSamples.percSignaturePOLE_MSI>20;
    elseif (strcmp(exclusionType, 'exclude_cll'))
        tableSamples.isExcluded(~ismember(tableSamples.project, {'MALY_DE', 'DLBC_US_LymphBNHL'})) = true;
    elseif (strcmp(exclusionType, 'exclude_lymphomas'))
        tableSamples.isExcluded(~ismember(tableSamples.project, {'CLLE_ES'})) = true;
    elseif (strcmp(exclusionType, 'excludeAsDig'))
        tableSamples.isExcluded(ismember(tableSamples.project, {'LINC_JP', 'LICA_FR', 'PAEN_AU', 'PAEN_IT'})) = true;
    else
        tableSamples.isExcluded = false(nSamples, 1);
    end
    tableSamples.isExcluded = tableSamples.isExcluded | ~tableSamples.isUsedOnePerDonor;
    if (sum(~tableSamples.isExcluded)==0), error('All samples excluded.'); end

    tableMutations.isExcluded = tableSamples.isExcluded(tableMutations.iSample);
    tableMutations.isHighCADD = tableMutations.CADD_PHRED >= minCADD_PHRED; % OLD: tableSNVs.isHighCADD = tableSNVs.iBinCADD >= minBinCADD;
    if (verbose)
        fprintf('\n\n\n%s: %d non-unique samples (multiple per donor)\n', tissueName, sum(~tableSamples.isUsedOnePerDonor));
        fprintf('%s %s: CODING DRIVER MUTATIONS: %d mutations (%d in enhancers).\n', tissueName, biosampleABC, sum(tableMutations.isPCAWGDriver), sum(tableMutations.isPCAWGDriver & tableMutations.iUniqueEnhancer>0));
        fprintf('%d included samples vs %d excluded samples (%.1f%%).\n', sum(~tableSamples.isExcluded), sum(tableSamples.isExcluded), 100*mean(tableSamples.isExcluded));
    end
    %%
    toc
    createDir(fileparts(saveFileData));
    save(saveFileData, 'tableSamples', 'tableMutations');
else
    fprintf('Loading data from %s...\n', saveFileData);
    load(saveFileData, 'tableSamples', 'tableMutations');
end