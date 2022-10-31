function yAltRelative = plotTFBS_logos(tableMutations_candidate, tableMotifs, iRow, isMOTIFG, geneName)
%lstCols = {'candidateGenes', 'MOTIFG_scoreDiff', 'MOTIFG_motifNamePrefix', 'tissuePrint', 'expressionMedianWT', 'expressionThisMut', 'gene', 'MOTIFG_motifName', 'MOTIFG_scoreAlt', 'MOTIFG_scoreRef', 'isHighCADD', 'VAF', 'qtlVAF'};
%tableMutations_candidate(iRow,lstCols)
%
if (isMOTIFG)
    % 19
    motifName = tableMutations_candidate.MOTIFG_motifName{iRow};
    motifStart = tableMutations_candidate.MOTIFG_motifStart(iRow);
    motifEnd = tableMutations_candidate.MOTIFG_motifEnd(iRow);
    positionMutation = tableMutations_candidate.MOTIFG_positionMutation(iRow);
    isMinusStrand = tableMutations_candidate.MOTIFG_isMinusStrand(iRow);
else
    % 49
    motifName = tableMutations_candidate.MOTIFBR_motifName{iRow};
    motifStart = tableMutations_candidate.MOTIFBR_motifStart(iRow);
    motifEnd = tableMutations_candidate.MOTIFBR_motifEnd(iRow);
    positionMutation = tableMutations_candidate.MOTIFBR_positionMutation(iRow);
    isMinusStrand = tableMutations_candidate.MOTIFBR_isMinusStrand(iRow);
end
%
context_ref = upper(tableMutations_candidate.context50bp{iRow});
lengthContext = length(context_ref);
iPosMutation = ceil(lengthContext/2);
if (~strcmp(context_ref(iPosMutation), tableMutations_candidate.ref{iRow})), error('Ref bases do not match.'); end
context_alt = context_ref; context_alt(iPosMutation) = tableMutations_candidate.alt{iRow};
x1 = motifStart-tableMutations_candidate.pos1(iRow)+1;    % iPosMutation + x1 >= 1 --> x1 >= 1-iPosMutation
x2 = motifEnd-tableMutations_candidate.pos1(iRow);        % iPosMutation + x2 <= lengthContext --> x2 <= lengthContext-iPosMutation = iPosMutation-1
if ((x1 < 1-iPosMutation) || (x2 > iPosMutation-1))
    error('Motif outside our context');
end
% context_ref(iPosMutation + (x1:x2))
sequence_ref = context_ref(iPosMutation + (x1:x2));
sequence_alt = context_alt(iPosMutation + (x1:x2));
% positionMutation = x2 + 1;
strandText = '+';
if (isMinusStrand)
    sequence_ref = seqrcomplement(sequence_ref);
    sequence_alt = seqrcomplement(sequence_alt);
    strandText = '-';
end
%

% strsplit(tableMotifs.motifName{1877}, '→')

iMotif = find(contains(tableMotifs.motifName, ['>', motifName]));
if (isempty(iMotif))
    error('Motif not found.')
elseif (length(iMotif)>1)
    warning('Multiple motif matches')
    iMotif = iMotif(1);
end

% iMotif = find(strcmp(cellfun(@(x) x(2:length(motifName)), tableMotifs.motifName), ['>', motifName]));
% tmp2 = cellfun(@(x) x(2:min([end, length(motifName)+1])), tableMotifs.motifName, 'UniformOutput', false);
% iMotif = find(strcmp(tmp2, motifName));
% if (length(iMotif)~=1)
%     error('Motif not found.')
% end

yAltRelative = plotMotif(tableMotifs.motifMatrix{iMotif}', 10, 2, sequence_ref, sequence_alt, motifName, positionMutation);
title(sprintf('{\\it%s} in %s\n%s:%s %s > %s (%s)', geneName, tableMutations_candidate.tissuePrint{iRow}, ...
    tableMutations_candidate.chr{iRow}, num2sepNumStr(tableMutations_candidate.pos1(iRow)), sequence_ref(positionMutation), sequence_alt(positionMutation), strandText)); % , tableMutations_candidate.ref{iRow}, tableMutations_candidate.alt{iRow}
%         title(sprintf('%s %s\n%s:%s %s>%s', tableMutations_candidate.candidateGenes{iRow}, tableMutations_candidate.tissuePrint{iRow}, ...
%             tableMutations_candidate.chr{iRow}, num2sepNumStr(tableMutations_candidate.pos1(iRow)), tableMutations_candidate.ref{iRow}, tableMutations_candidate.alt{iRow}));
end