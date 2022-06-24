function saveForOneGeneVisualisation(tissueName, biosampleABC, geneName, gene_pM, gene_pE, gene_qCombined, tableSamples, matCNV_genesSamples, matExpressionGenesSamples, matGenesSamplesNMut_SNVs_highCADD, ...
    tableMutations, matMutationsEnhancers, iGene, tableGencodeGenes, tableGenesNasserExpressed, matUniqueEnhancersGenes, tableUniqueEnhancers, tableUE_annotations_hyperUE, tableTrinucleotides)

CNVperSample=matCNV_genesSamples(iGene, :)';
expressionPerSample = matExpressionGenesSamples(iGene, :)';
isMutatedSample = matGenesSamplesNMut_SNVs_highCADD(iGene, :)'>0;

if (isempty(expressionPerSample))
    fprintf('No samples with expression for gene %s.\n', geneName);
    return
end

nSamples = size(tableSamples, 1);


sampleGroup = isMutatedSample + 1;
sampleGroup(tableSamples.isExcluded) = NaN;
expressionPerSample(tableSamples.isExcluded) = NaN;

sampleGroupInclWoExpression = sampleGroup;
sampleGroup(isnan(expressionPerSample)) = NaN;

if (sum(sampleGroup==2 & ~isnan(expressionPerSample))==0)
    fprintf('No mutated samples with expression for gene %s (%d mutated samples without expression).\n', geneName, sum(sampleGroup==2));
    return
end

%%

gene_pos0 = tableGencodeGenes.pos0(tableGenesNasserExpressed.iGencode(iGene));
gene_pos1 = tableGencodeGenes.pos1(tableGenesNasserExpressed.iGencode(iGene));
gene_TSS = tableGenesNasserExpressed.TSS(iGene);
gene_strand = tableGenesNasserExpressed.strand_GENCODE(iGene);

%%
% function ploGene_genomicView(imagesPath, tissueName, biosampleABC, geneName, pM, pE, FDR, tableSamples, isMutatedSample, expressionPerSample, CNVperSample, ...
%     iGene, tableGenesNasserExpressed, matUniqueEnhancersGenes, tableMutations_candidateOneTissue, ...
%     tableUniqueEnhancers, tableGencodeGenes, tableTrinucleotides, tableUE_annotations_hyperUE)
%
% nSamples = size(tableSamples, 1);


% imagesPathCurrent = [imagesPath, '/genomicView/'];

lstUniqueEnhancers = find(matUniqueEnhancersGenes(:,iGene)); %unique(tableEnhancers.iUniqueEnhancer(tableEnhancers.iGene == iGeneNasser));
gene_nUEs = length(lstUniqueEnhancers);

tableUniqueEnhancers_oneGene = tableUniqueEnhancers(lstUniqueEnhancers,:);
tableUE_annotations_hyperUE_oneGene = tableUE_annotations_hyperUE(lstUniqueEnhancers,:);


tableUE_annotations_hyperUE_oneGene.foldChangeScoreM = 2.^(tableUE_annotations_hyperUE_oneGene.eM_fullModel_SNVs_highCADD);
tableUE_annotations_hyperUE_oneGene.foldChangeScoreM(isinf(tableUE_annotations_hyperUE_oneGene.eM_fullModel_SNVs_highCADD)) = 0;

isMutPerEnhancer = full(matMutationsEnhancers(:,lstUniqueEnhancers)==1);
isMutOfThisGene = sum(isMutPerEnhancer, 2)>0;

isMutOK = isMutOfThisGene & tableMutations.isHighCADD & ~tableMutations.isExcluded;

tableMutationsThisGene = tableMutations(isMutOK,:);
isMutPerEnhancer = isMutPerEnhancer(isMutOK,:);

if (min(tableMutationsThisGene.iUniqueEnhancer)==0), error('Should be positive'); end
nPatterns = size(tableTrinucleotides, 1);
tableMutationsThisGene.expression = expressionPerSample(tableMutationsThisGene.iSample);
tableMutationsThisGene.yValues = log2(1+tableMutationsThisGene.expression);
tableMutationsThisGene.iPattern(tableMutationsThisGene.isIndel) = nPatterns;
tableMutationsThisGene.patternName = tableTrinucleotides.patternName(tableMutationsThisGene.iPattern);

save(['save/oneGene_', tissueName, '_', biosampleABC, '_', geneName], 'gene_pM', 'gene_pE', 'gene_qCombined', 'isMutatedSample', 'expressionPerSample', 'CNVperSample', 'nSamples', 'sampleGroup', 'sampleGroupInclWoExpression', ...
    'gene_pos0', 'gene_pos1', 'gene_TSS', 'gene_strand', 'gene_nUEs', 'tableMutationsThisGene', 'isMutPerEnhancer', 'tableUniqueEnhancers_oneGene', 'tableUE_annotations_hyperUE_oneGene');
%%
% fig = createMaximisedFigure(4);
% 
% yValues = log2(1+expressionPerSample);
% maxVal = max(yValues)*1.1;
% minVal = min(yValues)*0.9;
% cmap = lines(gene_nUEs);
% 
% x1 = sort([gene_pos0, gene_pos1]);
% 
% lstMarkers = {'o', 's', 'd', 'h'};
% 
% hLeg = []; hLegText = {'TSS'};
% hold on;
% hLeg(1) = plot(gene_TSS*[1,1], [minVal, maxVal], '-', 'LineWidth', 2, 'Color', 'y');
% if (strcmp(gene_strand, '+'))
%     text(gene_TSS, 0.95*maxVal, sprintf(' >> %s >>  >>  >>  >> ', geneName), 'Color', 'y', 'HorizontalAlignment', 'left');
% else
%     text(gene_TSS, 0.95*maxVal, sprintf(' <<  <<  <<  << %s << ', geneName), 'Color', 'y', 'HorizontalAlignment', 'right');
% end
% 
% y_b = minVal*[1,1];
% y_t = maxVal*[1,1];
% 
% for jUE = 1:gene_nUEs
%     colour = cmap(jUE, :);
%     %     iUE = lst_iUE(iColour);
%     %     isOK = tableMutationsThisGene.iUniqueEnhancer == iUE;
%     isOK = isMutPerEnhancer(:,jUE);
%     plot(tableMutationsThisGene.pos1(isOK), tableMutationsThisGene.yValues(isOK), lstMarkers{2}, 'MarkerFaceColor', colour, 'Color', (colour+1)/2);
% 
%     x2 = [tableUniqueEnhancers_oneGene.min_pos0(jUE), tableUniqueEnhancers_oneGene.max_pos1(jUE)];
%     h = patch([x2, fliplr(x2)], [y_b, fliplr(y_t)], colour, 'EdgeColor', colour);    set(h, 'FaceAlpha', 0.1);
% 
%     %         th = linspace( pi/2, -pi/2, 100);
%     %         R = 1;  %or whatever radius you want
%     %         x = R*cos(th) + 5;
%     %         y = R*sin(th) + 4;
%     %         plot(x,y);
% 
%     % BETTER BUT NOT COMPUTED ATM: label = {sprintf('%s, p=%s, mut=%d', tableUniqueEnhancers.name{iUE}, getPValueAsText(tableUniqueEnhancers.pBinomTestRight_SNVs_highCADD(iUE)), sum(isOK))};
%     label = {sprintf('%s, mut=%d, FC:%.1fx', tableUniqueEnhancers_oneGene.name{jUE}, sum(isOK), tableUE_annotations_hyperUE_oneGene.foldChangeScoreM(jUE))};
%     hLeg = [hLeg, h]; hLegText = [hLegText, label];
% end
% isOK = tableMutationsThisGene.isHighCADD;
% if (sum(isOK)>0)
%     h = plot(tableMutationsThisGene.pos1(isOK), tableMutationsThisGene.yValues(isOK), lstMarkers{4}, 'MarkerSize', 10, 'Color', .2*[1,1,1]);
%     hLeg = [hLeg, h]; hLegText = [hLegText, {'high CADD'}];
% end
% isOK = tableMutationsThisGene.isIndel;
% if (sum(isOK)>0)
%     h = plot(tableMutationsThisGene.pos1(isOK), tableMutationsThisGene.yValues(isOK), lstMarkers{3}, 'MarkerSize', 10, 'Color', .2*[1,1,1]);
%     for iMut = find(isOK)'
%         plot([tableMutationsThisGene.pos0(iMut), tableMutationsThisGene.pos1(iMut)], [tableMutationsThisGene.yValues(iMut), tableMutationsThisGene.yValues(iMut)], '-', 'MarkerSize', 10, 'Color', .2*[1,1,1]);
%     end
%     hLeg = [hLeg, h]; hLegText = [hLegText, {'indel'}];
% end
% 
% yLimVal = [minVal, maxVal]; clear h
% for iMut = find(isnan(tableMutationsThisGene.yValues))'
%     % OLD: [~, jUE] = ismember(tableMutationsThisGene.iUniqueEnhancer(iMut), lst_iUE); 
%     jUE = find(isMutPerEnhancer(iMut,:));
%     if (length(jUE)>1)
%         error('One mutation in multiple enhancers!');
%     end
%     %if (tmp2_mutations.isIndel(iMut)), marker = lstMarkers{3};
%     %elseif (tmp2_mutations.isHighCADD(iMut)), marker = lstMarkers{4}; else, marker = lstMarkers{2}; end [marker, ':']
%     h = plot(tableMutationsThisGene.pos1(iMut)*[1,1], yLimVal, ':', 'Color', cmap(jUE,:));
% end
% if (exist('h', 'var'))
%     hLeg = [hLeg, h]; hLegText = [hLegText, {'sample without expression'}];
% end
% h = plot(get(gca, 'XLim'), median(yValues(sampleGroup == 1), 'omitnan')*[1,1], '--r');
% hLeg = [hLeg, h]; hLegText = [hLegText, {'median WT expression'}];
% 
% %OLD legend([{'TSS'}; tableUniqueEnhancers.name(lst_iUE); {'high CADD'; 'known driver'; 'median WT'}], 'Location', 'EastOutside');
% 
% 
% text(tableMutationsThisGene.pos1, tableMutationsThisGene.yValues, tableMutationsThisGene.patternName);
% ylabel('Expression'); xlabel(sprintf('Position chr%d', tableMutationsThisGene.chrNumeric(1))); ylim(yLimVal);
% 
% 
% scaleLength = 20e3;
% %xStart = scaleLength*round(((minPos+maxPos)/2 - scaleLength/2)/scaleLength); yMiddle = yLimVal(2)*0.8; plot(xStart + [0,scaleLength], yMiddle*[1,1], '-k', 'LineWidth', 2);
% xStart = gene_TSS; yMiddle = yLimVal(2)*0.8; plot(xStart + [0,scaleLength], yMiddle*[1,1], '-k', 'LineWidth', 2);
% text(xStart, yMiddle*1.01, {sprintf('%d kbp', scaleLength/1e3), ''}, 'HorizontalAlignment', 'center');
% 
% xTickVal = get(gca, 'XTick');
% set(gca, 'XTickLabel', strcat(num2sepNumStr(round(xTickVal/1e3)), 'k'));
% 
% 
% xLimVal = get(gca, 'XLim');
% 
% h = patch([x1, fliplr(x1)], [0.90*maxVal*[1,1], fliplr(0.92*maxVal*[1,1])], [1,1,0], 'EdgeColor', 'none');    set(h, 'FaceAlpha', 0.5);   drawnow;
% %text(mean(x1), yLimVal(2)*.99, getPValueAsText(tableUniqueEnhancers.min_pBinomTestRight(iUE)), 'HorizontalAlignment', 'center', 'Rotation', 90);
% xlim(xLimVal);
% legend(hLeg, hLegText, 'Location', 'EastOutside');
% mySaveAs(fig, imagesPath, [tissueName, '_', geneName,biosampleABC,'.png']);


%% OLD
%         title(sprintf('%s: {\\bf%s} %s %s\n SNVs: %d (%d high-CADD) samples (%d mut total) | indels: %d samples (%d mut total) | {\\itpMut=%s} | {\\itpExpAdj=%s}', tissueName, geneNameSymbol, geneDriverText, strrep(biosampleABC, '_', ''),... % \nchr%d:%s-%s
%             length(unique(tmp2_mutations.iSample(~tmp2_mutations.isIndel))), length(unique(tmp2_mutations.iSample(tmp2_mutations.isHighCADD))), length(tmp2_mutations.iSample), ...
%             length(unique(tmp2_mutations.iSample(tmp2_mutations.isIndel))), length(unique(tmp2_mutations.iSample(tmp2_mutations.isIndel))), ...
%             getPValueAsText(tableGenesNasser.(xSignType)(iGeneNasser)), getPValueAsText(pValueExpressionAdjusted))); % tableGenesNasser.expression_pValue_min(iGeneNasser)