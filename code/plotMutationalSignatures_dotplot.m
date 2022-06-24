function plotMutationalSignatures_dotplot(tableMutations_candidate, tableTissues, sProperties)

nTissues = size(tableTissues, 1);

tableTrinucleotides = readtable(sProperties.TABLE_TRINUCLEOTIDES); % 'data/tableTriNucl96.txt'
nPatterns = size(tableTrinucleotides, 1);
tableSignaturesCOSMIC = readtable(sProperties.COSMIC_SIGNATURES_SBS); % 'data/COSMIC_v3.2_SBS_GRCh37.txt'
tableSignaturesCOSMIC.refBase = cellfun(@(x) x(3), tableSignaturesCOSMIC.Type, 'UniformOutput', false);             % C from 'G[C>G]G'
tableSignaturesCOSMIC.altBase = cellfun(@(x) x(5), tableSignaturesCOSMIC.Type, 'UniformOutput', false);             % G from 'G[C>G]G'
tableSignaturesCOSMIC.trinucleotide = cellfun(@(x) x([1,3,7]), tableSignaturesCOSMIC.Type, 'UniformOutput', false); % GCG from 'G[C>G]G'
tableSignaturesCOSMIC = sortrows(tableSignaturesCOSMIC,{'refBase', 'altBase', 'trinucleotide'},'ascend');
tableSignaturesCOSMIC.patternName = strcat(tableSignaturesCOSMIC.trinucleotide, '>', tableSignaturesCOSMIC.altBase);
if (~isequal(tableSignaturesCOSMIC.patternName, tableTrinucleotides.patternName)), error('WRONG ORDER'); end
tableSignaturesCOSMIC = tableSignaturesCOSMIC(:,2:end-4);
lstSignatures = tableSignaturesCOSMIC.Properties.VariableNames;
nSignatures = length(lstSignatures);
%%

matSimilarityTissuesSignatures = zeros(nTissues, nSignatures);

for iTissue = 1:nTissues
    isOK = tableMutations_candidate.iTissue == iTissue & ~tableMutations_candidate.isExcluded & tableMutations_candidate.isHighCADD;
    catalogue = histcounts(tableMutations_candidate.iPattern(isOK), 1:nPatterns+1);

    matSimilarityTissuesSignatures(iTissue,:) = computeSimilarityWithSignatures(tableSignaturesCOSMIC, catalogue');    
end
%%
% for iTissue = 1:nTissues
%     subplot(nR,nC,iTissue);
% 
%     isOK = tableMutations_candidate.iTissue == iTissue & ~tableMutations_candidate.isExcluded & tableMutations_candidate.isHighCADD;
%     vector = histcounts(tableMutations_candidate.iPattern(isOK), 1:nPatterns+1);
% 
%     plotMutationalSignature(vector, tableTissues.tissuePrint{iTissue}, imagesPath, ['Profile_', tableTissues.tissuePrint{iTissue}], tableTrinucleotides);
% end
%%

isSignAbove50 = max(matSimilarityTissuesSignatures, [], 1)>0.5;
matSimilarityTissuesSignatures_above50 = matSimilarityTissuesSignatures(:,isSignAbove50);
nSignAbove50 = sum(isSignAbove50);

minValue = min(matSimilarityTissuesSignatures_above50(:));
maxValue = max(matSimilarityTissuesSignatures_above50(:));

% fig = createMaximisedFigure(7, [0 0 35 7]); axes('Position', [0.07, 0.2, 0.9, 0.75]); 

hold on; cmap = flipud(lbmap(101, 'RedBlue')); fontSizeLarge = 12; fontSizeSmall = 8;

iLeg = 1; hLeg = NaN*ones(10,1);  legValues = cell(10,1); 
for value = linspace(maxValue, minValue, 10)
    valueColour = 1+round(100*(value-minValue)/(maxValue-minValue));
    colour = cmap(valueColour,:);
    hLeg(iLeg) = plot(-1, -1, 'o', 'MarkerSize', 1+floor(valueColour/5), 'Color', colour, 'MarkerFaceColor', (1+colour)/2);
    legValues{iLeg} = sprintf('%.1f', value);
    iLeg = iLeg + 1;
end

for iTissue = 1:nTissues
    for jSignature = 1:nSignAbove50
        value = matSimilarityTissuesSignatures_above50(iTissue, jSignature);
        valueColour = 1+round(100*(value-minValue)/(maxValue-minValue));
        colour = cmap(valueColour,:);
        plot(jSignature, iTissue, 'o', 'MarkerSize', 1+floor(valueColour/5), 'Color', colour, 'MarkerFaceColor', (1+colour)/2); 
        if (value > 0.75)
            text(jSignature, iTissue, sprintf('%.1f', value), 'Rotation', 0, 'HorizontalAlignment','center', 'FontSize', fontSizeSmall);
        end
    end
    drawnow;
end

legend(hLeg, legValues, 'Location', 'EastOutside'); ylim(0.5+[0,nTissues]);

xlabel('SBS mutational signatures'); xlim(0.5+[0,nSignAbove50]);
set(gca, 'YDir', 'reverse', 'YTick', 1:nTissues, 'YTickLabel', tableTissues.tissuePrint, 'XTick', 1:nSignAbove50, 'XTickLabel', strrep(lstSignatures(isSignAbove50), 'SBS', ''), 'XTickLabelRotation', 0, 'TickLength', [0 0], 'FontSize', fontSizeLarge);
% ax = gca;
% ax.XAxis.FontSize = fontSizeLarge;
% ax.YAxis.FontSize = fontSizeLarge;
% mySaveAs(fig, imagesPath, 'Fig3_signatures', false, false);
% %%
% vector = tableSignaturesCOSMIC.SBS84;
% textToShow = 'SBS84';
% saveNameFile = 'Fig3_signatureSBS84';
% plotMutationalSignature(vector, textToShow, imagesPath, saveNameFile, tableTrinucleotides)
% %%
% lstGenesTopHits = unique(tableMutations_candidate.candidateGenes(tableMutations_candidate.iTissue>1 & ~tableMutations_candidate.isExcluded & tableMutations_candidate.isHighCADD));
% nGenesTopHits = length(lstGenesTopHits);
% matSimilarityGenesSignatures = zeros(nGenesTopHits, nSignatures);
% 
% fig = createMaximisedFigure(4);
% nR = round(sqrt(nGenesTopHits)); nC = ceil(nGenesTopHits/nR); 
% 
% for iGeneTH = 1:nGenesTopHits
%     subplot(nR,nC,iGeneTH);
% 
% 
%     isOK = strcmp(tableMutations_candidate.candidateGenes, lstGenesTopHits{iGeneTH}) & ~tableMutations_candidate.isExcluded & tableMutations_candidate.isHighCADD;
% 
%     if (sum(isOK)>0)
%         catalogue = histcounts(tableMutations_candidate.iPattern(isOK), 1:nPatterns+1);
% 
%         matSimilarityGenesSignatures(iGeneTH,:) = computeSimilarityWithSignatures(tableSignaturesCOSMIC, catalogue');
% 
%         xValues = 1:nPatterns;
%         bar(xValues, catalogue);
% %         isOK = catalogue>0;
% %         text(xValues(isOK), catalogue(isOK), tableTrinucleotides.patternName(find(isOK)),'Rotation',90);
% %         ylim([0, ceil(max(catalogue(isOK))*1.2)]);
%         [maxVal, iMaxVal] = max(matSimilarityGenesSignatures(iGeneTH,:));
%         colour = 'k';
%         if (maxVal>.7)
%             colour = 'r';
%         end
%         set(gca, 'XTick', [], 'YTick', []);
%         title(sprintf('%s %.1f\n%s', lstGenesTopHits{iGeneTH}, maxVal, lstSignatures{iMaxVal}), 'Color', colour); drawnow;
%     end
% end
% mySaveAs(fig, imagesPath, 'Fig3_signaturesGenesSolid', false, false);
% %%
% % fig= createMaximisedFigure(1);
% % bar(matSimilarityGenesSignatures(:,68));
% % %%
% 
% %%
% lstGenesTopHits2 = lstGenesTopHits(matSimilarityGenesSignatures(:,68)<.3);
% nGenesTopHits2 = length(lstGenesTopHits2);
% matSimilarityGenesSignatures = zeros(nGenesTopHits2, nSignatures);
% 
% fig = createMaximisedFigure(5);
% nR = round(sqrt(nGenesTopHits2)); nC = ceil(nGenesTopHits2/nR); 
% 
% for iGeneTH = 1:nGenesTopHits2
%     subplot(nR,nC,iGeneTH);
% 
% 
%     isOK = strcmp(tableMutations_candidate.candidateGenes, lstGenesTopHits2{iGeneTH}) & ~tableMutations_candidate.isExcluded & tableMutations_candidate.isHighCADD;
% 
%     if (sum(isOK)>0)
%         catalogue = histcounts(tableMutations_candidate.iPattern(isOK), 1:nPatterns+1);
% 
%         matSimilarityGenesSignatures(iGeneTH,:) = computeSimilarityWithSignatures(tableSignaturesCOSMIC, catalogue');
% 
%         xValues = 1:nPatterns;
%         bar(xValues, catalogue);
% %         isOK = catalogue>0;
% %         text(xValues(isOK), catalogue(isOK), tableTrinucleotides.patternName(find(isOK)),'Rotation',90);
% %         ylim([0, ceil(max(catalogue(isOK))*1.2)]);
%         [maxVal, iMaxVal] = max(matSimilarityGenesSignatures(iGeneTH,:));
%         colour = 'k';
%         if (maxVal>.7)
%             colour = 'r';
%         end
%         set(gca, 'XTick', [], 'YTick', []);
%         title(sprintf('%s %.1f\n%s', lstGenesTopHits2{iGeneTH}, maxVal, lstSignatures{iMaxVal}), 'Color', colour); drawnow;
%     end
% end
% mySaveAs(fig, imagesPath, 'Fig3_signaturesGenesBlood_noSBS84', false, false);
% %% OLD: HEATMAP
% 
% %%
% % fig = createMaximisedFigure(6, [0 0 35 15]); cmap = flipud(lbmap(100, 'RedBlue')); colormap(cmap); fontSizeLarge = 20; fontSizeSmall = 10;
% % imagesc(matSimilarityTissuesSignatures); hB = colorbar; hB.Label.String = 'Cosine similarity'; hB.Label.FontSize = fontSizeLarge;
% % for iTissue = 1:nTissues
% %     for iSignature = 1:nSignatures
% %         if (matSimilarityTissuesSignatures(iTissue, iSignature) > 0.7)
% %             text(iSignature, iTissue, sprintf('%.2f', matSimilarityTissuesSignatures(iTissue, iSignature)), 'Rotation', 90, 'HorizontalAlignment','center', 'FontSize', fontSizeSmall);
% %         end
% %     end
% % end
% % caxis([0,1]);
% % hAxes = set(gca, 'YTick', 1:nTissues, 'YTickLabel', tableTissues.tissuePrint, 'XTick', 1:nSignatures, 'XTickLabel', lstSignatures, 'XTickLabelRotation', 90, 'TickLength', [0 0], 'FontSize', fontSizeLarge);
% % 
% % ax = gca;
% % ax.XAxis.FontSize = fontSizeSmall;
% % ax.YAxis.FontSize = fontSizeLarge;
% % 
% % mySaveAs(fig, imagesPath, 'Fig3_signatures', false, false);