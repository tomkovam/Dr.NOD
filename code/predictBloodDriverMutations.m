% function predictBloodDriverMutations(tableMutations_candidateOneTissue, geneName, tissueName, biosampleABC, exclusionType)
% NOT USED IN THE END

geneName = 'BCL2'; % BCL2 MYC

dataGene = load(['save/oneGene/oneGene_', tissueName, '_', biosampleABC, '_', geneName, '_', exclusionType], 'gene_pM', 'expressionPerSample', 'sampleGroup',  ...
    'gene_pos0', 'gene_pos1', 'gene_TSS', 'gene_strand', 'gene_nUEs', 'tableMutationsThisGene', 'isMutPerEnhancer', 'tableUniqueEnhancers_oneGene', 'tableUE_annotations_hyperUE_oneGene');

%%

%%
tmp = tableMutations_candidateOneTissue(contains(tableMutations_candidateOneTissue.candidateGenes, geneName) & ~tableMutations_candidateOneTissue.isExcluded,:);
tmp.expressionFC = tmp.expressionThisMut./tmp.expressionMedianWT; % , 'alt' , 'alt', 'CADD_PHRED', 'iPattern'
tmp2 = grpstats(tmp(:,{'chrNumeric', 'pos1', 'CADD_PHRED', 'iPattern', 'VAF', 'qtlVAF', 'FunSeq2_isAnnotated', 'FunSeq2_isMOTIFBR', 'FunSeq2_isMOTIFG', 'FunSeq2_gerp', 'FunSeq2_noncodingScore', 'isNearTSS_250bp', 'expressionMedianWT', 'expressionThisMut', 'expressionFC'}), ...
    {'chrNumeric', 'pos1'}, 'median');
tmp2.sortBy = tmp2.GroupCount.*tmp2.median_expressionFC;
tmp2 = sortrows(tmp2,'GroupCount','descend');
tmp3 = tmp2((tmp2.median_FunSeq2_isMOTIFBR | tmp2.median_FunSeq2_isMOTIFG) & tmp2.median_CADD_PHRED>=10 & tmp2.median_expressionFC > 2,:);
tmp3 = tmp2;
tmp3.printName = strrep(strcat(num2sepNumStr(tmp3.pos1), '_', num2str(tmp3.GroupCount), 's_', num2str(tmp3.median_CADD_PHRED, '%.0f'), 'CADD_', num2str(100*tmp3.median_VAF, '%.0f%%')), ' ', '');
%
nPositions = size(tmp3, 1);
nSamples = length(dataGene.expressionPerSample);
matPositionsSamples = false(nPositions, nSamples);
for iPosition = 1:nPositions
    matPositionsSamples(iPosition, tmp.iSample(tmp.pos1 == tmp3.pos1(iPosition))) = true;
end
%
expressionPerSample = dataGene.expressionPerSample;
isWTSample = dataGene.sampleGroup == 1;
isMutatedSample = dataGene.sampleGroup == 2;
expressionWT_median = median(expressionPerSample(isWTSample));
expressionWT_qtl75 = quantile(expressionPerSample(isWTSample), .75);
expressionWT_qtl25 = quantile(expressionPerSample(isWTSample), .25);
isExpressionAboveWT_median = expressionPerSample > expressionWT_median;
isExpressionAboveWT_qtl75 = expressionPerSample > expressionWT_qtl75;
isExpressionAboveWT_qtl25 = expressionPerSample > expressionWT_qtl25;
isExpressionAboveWT_median_2x = expressionPerSample > 2*expressionWT_median;
isExpressionAboveWT_qtl75_3x = expressionPerSample > 3*expressionWT_qtl75;
iSample = (1:nSamples)';
tableSamplesGene = table(iSample, expressionPerSample, isMutatedSample, isWTSample, isExpressionAboveWT_median, isExpressionAboveWT_qtl25, isExpressionAboveWT_qtl75, isExpressionAboveWT_median_2x, isExpressionAboveWT_qtl75_3x);
%
% fig = createMaximisedFigure(3);
% [~, permSamples] = sort(dataGene.expressionPerSample, 'descend');
% imagesc(matPositionsSamples(:,isSampleAboveWT)); colorbar;
tableSamplesGene.minPos = NaN*ones(nSamples, 1);
for iSample = find(isMutatedSample)'
    minPos = find(matPositionsSamples(:,iSample), 1, 'first');
    if (~isempty(minPos))
        tableSamplesGene.minPos(iSample) = minPos;
    end
end
tmp3.isSelected = false(nPositions, 1);
tmp3.isSelected(tableSamplesGene.minPos(tableSamplesGene.isMutatedSample & tableSamplesGene.isExpressionAboveWT_qtl75_3x & ~isnan(tableSamplesGene.minPos))) = true;
tmp3.isSelected = tmp3.pos1 >= 60986417 & tmp3.pos1 <= 60986433;
tmp3.isSelected = ismember(tmp3.pos1, [60986420, 60986421, 60986422, 60986423, 60986425, 60986426, 60986427, 60986430, 60986433]);
%
fig = createMaximisedFigure(4); hold on;

[~, permSamples] = sort(tableSamplesGene.expressionPerSample);
tableSamplesGene.indexSample(permSamples) = (1:nSamples)';

% xValues = (1:nSamples)';
xValues = tableSamplesGene.indexSample;
yValues = tableSamplesGene.expressionPerSample;
xValues(isnan(yValues)) = NaN;
plot(xValues, yValues, 'ok');
% text(xValues, yValues, num2str(tableSamplesGene.minPos));
isOK = dataGene.sampleGroup == 1;
plot(xValues(isOK), yValues(isOK), 'ok', 'MarkerFaceColor', .5*[1,1,1]);
lstPos = find(tmp3.isSelected);
nPos = length(lstPos); cmap = lines(nPos); hLeg = NaN*ones(nPos, 1);
for jPos = 1:nPos
    isOK = tmp.iSample(tmp.pos1 == tmp3.pos1(lstPos(jPos)));
    hLeg(jPos) = plot(xValues(isOK), yValues(isOK), 'o', 'MarkerSize', 3*jPos+2, 'Color', cmap(jPos,:), 'LineWidth', 2);
end
legend(hLeg, tmp3.printName(lstPos), 'Interpreter', 'none', 'Location', 'NorthWest');
xlim([1, max(xValues)]);
%%
tableSamplesGene.hasNR3C1_break = false(nSamples, 1);
tableSamplesGene.hasNR3C1_break(tmp.iSample(ismember(tmp.pos1, tmp3.pos1(tmp3.isSelected)))) = true;
fig = createMaximisedFigure(6);
yValues = tableSamplesGene.expressionPerSample;
gValues = tableSamplesGene.isMutatedSample + 1;
gValues(tableSamplesGene.hasNR3C1_break) = 3;
boxchart(gValues, yValues);
tableSamplesGene(gValues==2 & yValues>70,:)
sum(gValues==3 & yValues>expressionWT_qtl75)
sum(gValues==2 & yValues>expressionWT_qtl75)
sum(yValues>expressionWT_qtl75)
%%
tableFunSeq2_a = readtable('data/FunSeq2/annotatedPCAWG/intersectedB.PCAWG_DLBC_US.motif.FunSeq2.bed.txt');
tableFunSeq2_b = readtable('data/FunSeq2/annotatedPCAWG/intersectedB.PCAWG_MALY_DE.motif.FunSeq2.bed.txt');
tableFunSeq2 = [tableFunSeq2_a; tableFunSeq2_b];
%%
% tmpSample = tmp(tmp.iSample == 150,:);
tmpSample = tmp(ismember(tmp.iSample, tableSamplesGene.iSample(gValues == 2 & yValues > 2*expressionWT_qtl75)),:);
tmpFS = unique(tableFunSeq2(ismember(tableFunSeq2.Var3, tmpSample.pos1) & strcmp(tableFunSeq2.Var1, sprintf('chr%d', tmpSample.chrNumeric(1))),:));
tmpFS = sortrows(tmpFS,'Var3','descend');
tmpFS.isTFBS_event = ~strcmp(tmpFS.Var13, '.');
[isOK, index] = ismember(tmpFS.Var3, tmp3.pos1); mean(isOK)
tmpFS.index_tmp3 = index;
tmpFS.mutatedSamples = tmp3.GroupCount(index);
fig = createMaximisedFigure(7);
isOK = tmpFS.isTFBS_event;
plot(tmpFS.Var3(isOK), tmpFS.mutatedSamples(isOK), 'o');
xlim([60984500, 60988800]);
xlim([60988040, 60988060]);
xlim([60988220, 60988240]);
xlim([60988220, 60988250]);
tmpFS.pos1_suffix = tmpFS.Var3 - 60988000;
tmpFS_selection = tmpFS(tmpFS.Var3>=60988220 & tmpFS.Var3<=60988250, {'mutatedSamples', 'pos1_suffix', 'Var3', 'Var13'});
tmpFS_selection = sortrows(tmpFS_selection,'pos1_suffix','descend');
tmpFS_selection(contains(tmpFS_selection.Var13, '60988218#60988236'),:) % motif-break of FOXJ2, FOXO3, FOXF2, FOXI, FOXA, HDAC, EP300, NKX2
unique(tmpFS_selection.Var3(contains(tmpFS_selection.Var13, '60988218#60988236'),:))
tmpFS_selection(contains(tmpFS_selection.Var13, 'MOTIFG'),:)
% unique(tmpFS.Var3(tmpFS.mutatedSamples>5 & tmpFS.isTFBS_event,:))
%%
lstSecondCluster = [60988222    60988223    60988225    60988228    60988229    60988230    60988231    60988235];
tableSamplesGene.hasOther_break = false(nSamples, 1);
tableSamplesGene.hasOther_break(tmp.iSample(ismember(tmp.pos1, lstSecondCluster))) = true;
sum(tableSamplesGene.hasOther_break) % 16
sum(gValues==3 & yValues>expressionWT_qtl75 & tableSamplesGene.hasOther_break) % 7
sum(gValues==2 & yValues>expressionWT_qtl75 & tableSamplesGene.hasOther_break) % 9

min(yValues(gValues == 3))
%%
[~, ~, ~, ~, ~, ~, ~, ...
    ~, ~, ~, ~, ...
    ~, ~, ~, ~, matCNV_genesSamples, matSV_genesSamples_minDistance, matSV_genesSamples_nameSVs] = ...
    loadCancerData(runAgain, tissueName, biosampleABC, sProperties.enhancerAnalysis, false, true, tissueNameSV, sProperties);
%%

%%
% fig = createMaximisedFigure(2);
% plot(tmp3.GroupCount, tmp3.median_expressionFC, 'o')
%%
% fig = createMaximisedFigure(1);
% xValues = tmp.pos1;
% yValues = tmp.expressionFC;
% plot(xValues, yValues, 'o');