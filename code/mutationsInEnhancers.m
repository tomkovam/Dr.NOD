function [tableMutations, matMutationsEnhancers] = mutationsInEnhancers(tableMutations, tableUniqueEnhancers_regions, tableChrSizes)

% The input tableMutations needs to have these columns: chrNumeric, pos0, pos1
% For every mutation, we annotate it with iUniqueEnhancer.
% We go through the rows in tableUniqueEnhancers_regions and for each row we look at which iUniqueEnhancer it belongs to and we annotate all mutations in that region that they fall into iUniqueEnhancer.

nMutations = size(tableMutations, 1);
nUniqueEnhancers = max(tableUniqueEnhancers_regions.iUniqueEnhancer);
nChromosomes = size(tableChrSizes, 1);
%% For each enhancer, compute the number of variants inside and around the enhancer
t = tic;
fprintf('RUNNING mutationsInEnhancer...\n');

matMutationsEnhancers = sparse(nMutations, nUniqueEnhancers);
for iChr = 1:nChromosomes % By splitting it into chromosomes, we speed it up
    lstEnhancers = find(tableUniqueEnhancers_regions.chrNumeric == iChr)';
    tmp_tableMutations = tableMutations(tableMutations.chrNumeric == iChr,:);
    tmp_tableMutations.iUniqueEnhancer = zeros(size(tmp_tableMutations, 1), 1);
    tmp_matMutationsEnhancers = sparse(size(tmp_tableMutations, 1), nUniqueEnhancers);
    for iUniqueEnhancer_region = lstEnhancers
        iUniqueEnhancer = tableUniqueEnhancers_regions.iUniqueEnhancer(iUniqueEnhancer_region);
        isOK = tmp_tableMutations.pos1 > tableUniqueEnhancers_regions.pos0(iUniqueEnhancer_region) & tmp_tableMutations.pos0 < tableUniqueEnhancers_regions.pos1(iUniqueEnhancer_region);
        tmp_tableMutations.iUniqueEnhancer(isOK) = iUniqueEnhancer;
        tmp_matMutationsEnhancers(isOK, iUniqueEnhancer) = true;        
    end
    tableMutations.iUniqueEnhancer(tableMutations.chrNumeric == iChr) = tmp_tableMutations.iUniqueEnhancer;
    matMutationsEnhancers(tableMutations.chrNumeric == iChr, :) = tmp_matMutationsEnhancers;
end
toc(t)
%%