function [dataTFBS, tableMutations_candidate] = loadMotifs(tableMutations_candidate, tableTissues_data4, sProperties)

tableDataMotifs = readtable(sProperties.FUNSEQ2_ENCODE_MOTIFS, 'ReadVariableNames', false, 'NumHeaderLines', 0); % 'data/FunSeq2/ENCODE_motifs.txt' From http://compbio.mit.edu/encode-motifs/ | http://compbio.mit.edu/encode-motifs/motifs.txt
tableDataMotifs.Properties.VariableNames = {'consensus', 'A', 'C', 'G', 'T'};
tableDataMotifs.isHeaderLine = isnan(tableDataMotifs.A);
nMotifs = sum(tableDataMotifs.isHeaderLine);

indexHeader = find(tableDataMotifs.isHeaderLine);
indexHeader = [indexHeader; size(tableDataMotifs, 1)+1];

tableMotifs = table();
tableMotifs.motifName = cell(nMotifs, 1);
tableMotifs.motifMatrix = cell(nMotifs, 1);
tableMotifs.motifConsensus = cell(nMotifs, 1);
for iMotif = 1:nMotifs
    minRow = indexHeader(iMotif);
    maxRow = indexHeader(iMotif+1)-1;
    tableMotifs.motifName{iMotif} = tableDataMotifs.consensus{minRow};
    tableMotifs.motifMatrix{iMotif} = table2array(tableDataMotifs(minRow+1:maxRow, {'A', 'C', 'G', 'T'}));
    tableMotifs.motifConsensus{iMotif} = strjoin(tableDataMotifs.consensus(minRow+1:maxRow)', '');
end
%%
tableMotifs.motifPrefix = regexp(tableMotifs.motifName, '(?<=>)[^_]+(?=_)', 'match', 'once');
tableMotifPrefix = table();
tableMotifPrefix.motifPrefix = unique(tableMotifs.motifPrefix);
nMotifPrefix = length(tableMotifPrefix.motifPrefix);
%%
tableMutations_candidate.MOTIFG = arrayfun(@(x) regexp(x, 'MOTIFG=.+', 'match', 'once'), tableMutations_candidate.FunSeq2_motif_analysis, 'UniformOutput', true);
if (~isequal(~strcmp(tableMutations_candidate.MOTIFG, ''), tableMutations_candidate.isMOTIFG==1)), error('MOTIFBR assignment not correct'); end
if (sum(contains(tableMutations_candidate.MOTIFG, 'MOTIFBR'))>0), error('MOTIFBR in MOTIFBR'); end
%
nRows = size(tableMutations_candidate, 1);
tableMutations_candidate.MOTIFG_motifName = cell(nRows, 1); tableMutations_candidate.MOTIFG_motifName(:) = {''};
tableMutations_candidate.MOTIFG_motifNamePrefix = cell(nRows, 1); tableMutations_candidate.MOTIFG_motifNamePrefix(:) = {''};
tableMutations_candidate.MOTIFG_motifStart = NaN*ones(nRows, 1); 
tableMutations_candidate.MOTIFG_motifEnd = NaN*ones(nRows, 1); 
tableMutations_candidate.MOTIFG_isMinusStrand = false(nRows, 1); 
tableMutations_candidate.MOTIFG_positionMutation = NaN*ones(nRows, 1); 
tableMutations_candidate.MOTIFG_scoreAlt = NaN*ones(nRows, 1); 
tableMutations_candidate.MOTIFG_scoreRef = NaN*ones(nRows, 1); 
tableMutations_candidate.MOTIFG_scoreDiff = NaN*ones(nRows, 1); 
for iRow = 1:nRows
    if (~isempty(tableMutations_candidate.MOTIFG{iRow}))
        tmp2 = strsplit(tableMutations_candidate.MOTIFG{iRow}, {'=', '#'}); % , ','
        tableMutations_candidate.MOTIFG_motifName{iRow} = tmp2{2};
        tmp3 = regexp(tmp2{2}, '.*_', 'match', 'once');
        tableMutations_candidate.MOTIFG_motifNamePrefix{iRow} = tmp3(1:end-1);
        tableMutations_candidate.MOTIFG_motifStart(iRow) = str2double(tmp2{3});
        tableMutations_candidate.MOTIFG_motifEnd(iRow) = str2double(tmp2{4});
        tableMutations_candidate.MOTIFG_isMinusStrand(iRow) = strcmp(tmp2{5}, '-');
        tableMutations_candidate.MOTIFG_positionMutation(iRow) = str2double(tmp2{6});
        tableMutations_candidate.MOTIFG_scoreAlt(iRow) = str2double(tmp2{7});               % sequence score with alternative allele
        tmp3 = strsplit(tmp2{8}, ',');
        tableMutations_candidate.MOTIFG_scoreRef(iRow) = str2double(tmp3{1});               % sequence score with reference allele
    end
end
tableMutations_candidate.MOTIFG_scoreDiff = tableMutations_candidate.MOTIFG_scoreAlt - tableMutations_candidate.MOTIFG_scoreRef;
%
tableMutations_candidate.MOTIFBR = arrayfun(@(x) regexp(x, 'MOTIFBR=.+', 'match', 'once'), tableMutations_candidate.FunSeq2_motif_analysis, 'UniformOutput', true);
isOK = contains(tableMutations_candidate.FunSeq2_motif_analysis, 'MOTIFG');
tableMutations_candidate.MOTIFBR(isOK) = arrayfun(@(x) regexp(x, 'MOTIFBR=.+(?=MOTIFG)', 'match', 'once'), tableMutations_candidate.FunSeq2_motif_analysis(isOK), 'UniformOutput', true);
if (~isequal(~strcmp(tableMutations_candidate.MOTIFBR, ''), tableMutations_candidate.isMOTIFBR==1)), error('MOTIFBR assignment not correct'); end
if (sum(contains(tableMutations_candidate.MOTIFBR, 'MOTIFG'))>0), error('MOTIFG in MOTIFBR'); end
tableMutations_candidate.MOTIFBR_TFs = cell(nRows, 1); tableMutations_candidate.MOTIFBR_TFs(:) = {''};
tableMutations_candidate.MOTIFBR_motifName = cell(nRows, 1); tableMutations_candidate.MOTIFBR_motifName(:) = {''};
tableMutations_candidate.MOTIFBR_motifNamePrefix = cell(nRows, 1); tableMutations_candidate.MOTIFBR_motifNamePrefix(:) = {''};
tableMutations_candidate.MOTIFBR_motifStart = NaN*ones(nRows, 1); 
tableMutations_candidate.MOTIFBR_motifEnd = NaN*ones(nRows, 1); 
tableMutations_candidate.MOTIFBR_isMinusStrand = false(nRows, 1); 
tableMutations_candidate.MOTIFBR_positionMutation = NaN*ones(nRows, 1); 
tableMutations_candidate.MOTIFBR_scoreAlt = NaN*ones(nRows, 1); 
tableMutations_candidate.MOTIFBR_scoreRef = NaN*ones(nRows, 1); 
tableMutations_candidate.MOTIFBR_scoreDiff = NaN*ones(nRows, 1); 
for iRow = 1:nRows
    if (~isempty(tableMutations_candidate.MOTIFBR{iRow}))
        tmp2 = strsplit(tableMutations_candidate.MOTIFBR{iRow}, {'=', '#'}); % , ','
        tableMutations_candidate.MOTIFBR_TFs{iRow} = tmp2{2};
        tableMutations_candidate.MOTIFBR_motifName{iRow} = tmp2{3};
        tmp3 = regexp(tmp2{3}, '.*_', 'match', 'once');
        tableMutations_candidate.MOTIFBR_motifNamePrefix{iRow} = tmp3(1:end-1);
        tableMutations_candidate.MOTIFBR_motifStart(iRow) = str2double(tmp2{4});
        tableMutations_candidate.MOTIFBR_motifEnd(iRow) = str2double(tmp2{5});
        tableMutations_candidate.MOTIFBR_isMinusStrand(iRow) = strcmp(tmp2{6}, '-');
        tableMutations_candidate.MOTIFBR_positionMutation(iRow) = str2double(tmp2{7});
        tableMutations_candidate.MOTIFBR_scoreAlt(iRow) = str2double(tmp2{8});               % alternative allele frequency in PWM
        tmp3 = strsplit(tmp2{9}, ',');
        tableMutations_candidate.MOTIFBR_scoreRef(iRow) = str2double(tmp3{1});               % reference allele frequency in PWM
    end
end
tableMutations_candidate.MOTIFBR_scoreDiff = tableMutations_candidate.MOTIFBR_scoreRef - tableMutations_candidate.MOTIFBR_scoreAlt;
%%
isOK = ~tableMutations_candidate.isIndel & ~tableMutations_candidate.isExcluded & tableMutations_candidate.iTissue>1;% & (tableMutations_candidate.isMOTIFBR==1 | tableMutations_candidate.isMOTIFG==1);
tableMutations_candidate.isOK = isOK;
isMutMotifPrefix_MOTIFBR = false(size(tableMutations_candidate, 1), nMotifPrefix);
isMutMotifPrefix_MOTIFG = false(size(tableMutations_candidate, 1), nMotifPrefix);
for iMotifPrefix = 1:nMotifPrefix
    motifPrefix = tableMotifPrefix.motifPrefix{iMotifPrefix};
    isMutMotifPrefix_MOTIFBR(:, iMotifPrefix) = contains(tableMutations_candidate.MOTIFBR, ['=', motifPrefix, '_']) |...
        contains(tableMutations_candidate.MOTIFBR, [',', motifPrefix, '_'])|...
        contains(tableMutations_candidate.MOTIFBR, ['#', motifPrefix, '_']);
    isMutMotifPrefix_MOTIFG(:, iMotifPrefix) = contains(tableMutations_candidate.MOTIFG, ['=', motifPrefix, '_']) |...
        contains(tableMutations_candidate.MOTIFG, [',', motifPrefix, '_'])|...
        contains(tableMutations_candidate.MOTIFG, ['#', motifPrefix, '_']);
end
tableMutations_candidate.nMotifPrefix_MOTIFBR = sum(isMutMotifPrefix_MOTIFBR, 2);
tableMutations_candidate.nMotifPrefix_MOTIFG = sum(isMutMotifPrefix_MOTIFG, 2);
tableMotifPrefix.nMOTIFBR_UP = sum(isMutMotifPrefix_MOTIFBR(isOK & tableMutations_candidate.isCandidateDriverUP,:), 1)';
tableMotifPrefix.nMOTIFG_UP = sum(isMutMotifPrefix_MOTIFG(isOK & tableMutations_candidate.isCandidateDriverUP,:), 1)';
tableMotifPrefix.nMOTIFBR_DOWN = sum(isMutMotifPrefix_MOTIFBR(isOK & ~tableMutations_candidate.isCandidateDriverUP,:), 1)';
tableMotifPrefix.nMOTIFG_DOWN = sum(isMutMotifPrefix_MOTIFG(isOK & ~tableMutations_candidate.isCandidateDriverUP,:), 1)';
if (sum((tableMutations_candidate.isMOTIFBR==1 | tableMutations_candidate.isMOTIFG==1) & (tableMutations_candidate.nMotifPrefix_MOTIFBR+tableMutations_candidate.nMotifPrefix_MOTIFG)==0)>0)
    error('Some MOTIFBR/MOTIFBG do not have a motif prefix assigned')
end
if (sum(~(tableMutations_candidate.isMOTIFBR==1 | tableMutations_candidate.isMOTIFG==1) & (tableMutations_candidate.nMotifPrefix_MOTIFBR+tableMutations_candidate.nMotifPrefix_MOTIFG)>0)>0)
    error('Some nonMOTIFBR/MOTIFBG do have a motif prefix assigned')
end
%
isOK = tableMutations_candidate.isOK;
tmp = tableMutations_candidate(isOK,:);
tmp2_MOTIFBR = isMutMotifPrefix_MOTIFBR(isOK,:);
tmp2_MOTIFG = isMutMotifPrefix_MOTIFG(isOK,:);
nTissues = size(tableTissues_data4, 1);
matTissuesMotifPrefix_MOTIFBR_UP = zeros(nTissues, nMotifPrefix);
matTissuesMotifPrefix_MOTIFG_UP = zeros(nTissues, nMotifPrefix);
matTissuesMotifPrefix_MOTIFBR_DOWN = zeros(nTissues, nMotifPrefix);
matTissuesMotifPrefix_MOTIFG_DOWN = zeros(nTissues, nMotifPrefix);
for iTissue = 1:nTissues
    matTissuesMotifPrefix_MOTIFBR_UP(iTissue,:) = sum(tmp2_MOTIFBR(tmp.iTissue==iTissue & tmp.isMOTIFBR==1 & tmp.isCandidateDriverUP,:), 1);
    matTissuesMotifPrefix_MOTIFG_UP(iTissue,:) = sum(tmp2_MOTIFG(tmp.iTissue==iTissue & tmp.isMOTIFG==1 & tmp.isCandidateDriverUP,:), 1);
    matTissuesMotifPrefix_MOTIFBR_DOWN(iTissue,:) = sum(tmp2_MOTIFBR(tmp.iTissue==iTissue & tmp.isMOTIFBR==1 & ~tmp.isCandidateDriverUP,:), 1);
    matTissuesMotifPrefix_MOTIFG_DOWN(iTissue,:) = sum(tmp2_MOTIFG(tmp.iTissue==iTissue & tmp.isMOTIFG==1 & ~tmp.isCandidateDriverUP,:), 1);
end
%%
if (false)
    %%
    fig = createMaximisedFigure(3); nTissues = size(matTissuesMotifPrefix_MOTIFBR_UP, 1);
    subplot(2,2,1);
    isOK = tableMotifPrefix.nMOTIFBR_UP>3;
    imagesc(matTissuesMotifPrefix_MOTIFBR_UP(2:nTissues,isOK)); colorbar;
    set(gca, 'XTick', 1:sum(isOK), 'XTickLabel', tableMotifPrefix.motifPrefix(isOK), 'YTick', 1:nTissues-1, 'YTickLabel', tableTissues_data4.tissuePrint(2:nTissues));
    subplot(2,2,2);
    isOK = tableMotifPrefix.nMOTIFG_UP>1;
    imagesc(matTissuesMotifPrefix_MOTIFG_UP(2:nTissues,isOK)); colorbar;
    set(gca, 'XTick', 1:sum(isOK), 'XTickLabel', tableMotifPrefix.motifPrefix(isOK), 'YTick', 1:nTissues-1, 'YTickLabel', tableTissues_data4.tissuePrint(2:nTissues));
    subplot(2,2,3);
    isOK = tableMotifPrefix.nMOTIFBR_DOWN>0;
    imagesc(matTissuesMotifPrefix_MOTIFBR_DOWN(2:nTissues,isOK)); colorbar;
    set(gca, 'XTick', 1:sum(isOK), 'XTickLabel', tableMotifPrefix.motifPrefix(isOK), 'YTick', 1:nTissues-1, 'YTickLabel', tableTissues_data4.tissuePrint(2:nTissues));
    subplot(2,2,4);
    isOK = tableMotifPrefix.nMOTIFG_DOWN>0;
    imagesc(matTissuesMotifPrefix_MOTIFG_DOWN(2:nTissues,isOK)); colorbar;
    set(gca, 'XTick', 1:sum(isOK), 'XTickLabel', tableMotifPrefix.motifPrefix(isOK), 'YTick', 1:nTissues-1, 'YTickLabel', tableTissues_data4.tissuePrint(2:nTissues));
end
%% UNUSED (most are both positive and negative regulators)
if (false)
    tablePositiveRegulator = readtable(sProperties.GO_POSITIVE_REGULATION); % 'data/QuickGO/QuickGO-annotations-0045893-20220414.txt' positive regulation of transcription, DNA-templated
    tableNegativeRegulator = readtable(sProperties.GO_NEGATIVE_REGULATION); % 'data/QuickGO/QuickGO-annotations-0045892-20220414.txt' negative regulation of transcription, DNA-templated
    lstPositiveRegulator = unique(tablePositiveRegulator.SYMBOL);
    lstNegativeRegulator = unique(tableNegativeRegulator.SYMBOL);
    tableMotifPrefix.isPositiveRegulator = ismember(tableMotifPrefix.motifPrefix, lstPositiveRegulator);
    tableMotifPrefix.isNegativeRegulator = ismember(tableMotifPrefix.motifPrefix, lstNegativeRegulator);
    tableMotifPrefix.isPositiveRegulator_prefix = cellfun(@(x) max(contains(lstPositiveRegulator, x)), tableMotifPrefix.motifPrefix);
    tableMotifPrefix.isNegativeRegulator_prefix = cellfun(@(x) max(contains(lstNegativeRegulator, x)), tableMotifPrefix.motifPrefix);
end
%%
dataTFBS.tableMotifs = tableMotifs;
dataTFBS.tableMotifPrefix = tableMotifPrefix;
dataTFBS.isMutMotifPrefix_MOTIFBR = isMutMotifPrefix_MOTIFBR;
dataTFBS.isMutMotifPrefix_MOTIFG = isMutMotifPrefix_MOTIFG;
dataTFBS.matTissuesMotifPrefix_MOTIFBR_UP = matTissuesMotifPrefix_MOTIFBR_UP;
dataTFBS.matTissuesMotifPrefix_MOTIFG_UP = matTissuesMotifPrefix_MOTIFG_UP;
dataTFBS.matTissuesMotifPrefix_MOTIFBR_DOWN = matTissuesMotifPrefix_MOTIFBR_DOWN;
dataTFBS.matTissuesMotifPrefix_MOTIFG_DOWN = matTissuesMotifPrefix_MOTIFG_DOWN;
