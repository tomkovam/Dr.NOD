function [tableDataForBMM_averaged, tableDataForBMM_averaged_info] = slidingWindowEnhancers(tableUniqueEnhancers, nTheoreticalMutations_PHRED_geqCUTOFF, tableDataForBMM, nMutSamplesPerRow, maxSize)


% maxSize = 1300;
tableUniqueEnhancersCopy = tableUniqueEnhancers;

tableUniqueEnhancers.index = (1:size(tableUniqueEnhancers, 1))'; % To save the original order
[~, permRows] = sortrows(tableUniqueEnhancers(:,{'chrNumeric','min_pos0', 'max_pos1'}));
tableUniqueEnhancers = tableUniqueEnhancers(permRows,:);
tableDataForBMM = tableDataForBMM(permRows,:);
nMutSamplesPerRow = nMutSamplesPerRow(permRows,:);
nTheoreticalMutations_PHRED_geqCUTOFF = nTheoreticalMutations_PHRED_geqCUTOFF(permRows,:);
nRows = length(permRows);
%%
tableDataForBMM_averaged = tableDataForBMM;

matDataForBMM = table2array(tableDataForBMM_averaged);
matDataForBMM_averaged = NaN*matDataForBMM;

tableDataForBMM_averaged_info = table();
tableDataForBMM_averaged_info.min_pos0 = tableUniqueEnhancers.min_pos0;
tableDataForBMM_averaged_info.nAveraged = NaN*ones(nRows, 1);

tableUniqueEnhancers.nPositions_cumsum = cumsum(tableUniqueEnhancers.nPositions);
for iRow = 1:nRows
%     if (mod(iRow, 1000)==1)
%         fprintf('Processing row %d...\n', iRow);
%     end

    % By distance
    %     maxDistance = 1e7; % 1e7 | 2.9e5 --> median of 10,149 bp and explained variance=0.0798 | 1e6 --> median of 25,244 bp and explained variance=0.1639 | 1e7 --> median of 167,992 bp bp and explained variance=0.5919
    %     min_pos0 = tableUniqueEnhancers.min_pos0(iRow);
    %     max_pos1 = tableUniqueEnhancers.max_pos1(iRow);
    %     chrNumeric = tableUniqueEnhancers.chrNumeric(iRow);
    %     isOK = tableUniqueEnhancers.chrNumeric==chrNumeric & ...
    %         ((min_pos0 - tableUniqueEnhancers.max_pos1) < maxDistance) & ...
    %         ((tableUniqueEnhancers.min_pos0 - max_pos1) < maxDistance) & ...
    %         tableUniqueEnhancers.mean_blacklisted==0;
    
    % By size
    chrNumeric = tableUniqueEnhancers.chrNumeric(iRow);    
    if (iRow == 1)
        baseline = 0;
    else
        baseline = tableUniqueEnhancers.nPositions_cumsum(iRow-1);
    end
    tmp = tableUniqueEnhancers.nPositions_cumsum - baseline;
    isOK = tableUniqueEnhancers.chrNumeric==chrNumeric & tmp > 0 & tmp <= maxSize;
    %     isOK = false(nRows, 1); isOK(iRow) = true;
    %     if (iRow<nRows)
    %         isOK(iRow+1) = true;
    %     end
    %     if (iRow+1<nRows)
    %         isOK(iRow+2) = true;
    %     end
    
    %
    tableDataForBMM_averaged_info.nAveraged(iRow) = sum(isOK);
    nPositionsInEnhancersTotal = sum(tableDataForBMM.nPositionsInEnhancers(isOK));
    matDataForBMM_averaged(iRow,:) = mean(matDataForBMM(isOK,:),1,'omitnan'); % unweighted average
    %matDataForBMM_averaged(iRow,:) = sum(tableDataForBMM.nPositionsInEnhancers(isOK,:).*matDataForBMM(isOK,:),'omitnan')/nPositionsInEnhancersTotal; % weighted average (weighted by the size of each enhancer)
    tableDataForBMM_averaged_info.nPositionsInEnhancers(iRow) = nPositionsInEnhancersTotal;
    tableDataForBMM_averaged_info.nPositionsInEnhancersCADD(iRow) = sum(nTheoreticalMutations_PHRED_geqCUTOFF(isOK));
    tableDataForBMM_averaged_info.nTheoreticalMutations_PHRED_geqCUTOFF(iRow) = sum(nTheoreticalMutations_PHRED_geqCUTOFF(isOK));
    tableDataForBMM_averaged_info.mutationFrequency(iRow) = sum(nMutSamplesPerRow(isOK))/sum(nTheoreticalMutations_PHRED_geqCUTOFF(isOK));
end
tableDataForBMM_averaged_info.mutationFrequency(isinf(tableDataForBMM_averaged_info.mutationFrequency)) = NaN; % This can happen when nTheoreticalMutations_PHRED_geqCUTOFF=0 and nMutSamplesInEnhancersPerGene>0 because of indels (should not happen for highCADD SNVs)
lstPredictors = tableDataForBMM_averaged.Properties.VariableNames;
tableDataForBMM_averaged = array2table(matDataForBMM_averaged);
tableDataForBMM_averaged.Properties.VariableNames = lstPredictors;
tableDataForBMM_averaged.mutationFrequency = tableDataForBMM_averaged_info.mutationFrequency;

%%
[~, permRows] = sortrows(tableUniqueEnhancers(:,{'index'})); % To restore the original order
tableDataForBMM_averaged = tableDataForBMM_averaged(permRows,:);
tableDataForBMM_averaged_info = tableDataForBMM_averaged_info(permRows,:);
%%
% median(tableDataForBMM_averaged_info.nPositionsInEnhancers)
if (~isequal(tableUniqueEnhancersCopy.min_pos0, tableDataForBMM_averaged_info.min_pos0)), error('Rows sorted incorrectly.'); end
