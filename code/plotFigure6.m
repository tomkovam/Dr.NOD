function plotFigure6(tableMutationGenePairs, imagesPath, sColours, tableTissues_data1)

rng(1);
tmp1 = tableMutationGenePairs(tableMutationGenePairs.isHighCADD & ~tableMutationGenePairs.isExcluded, :);

nTissues = size(tableTissues_data1, 1);
fontSize = 12;

for iType = 1:2
    if (iType == 1)
        isCloserToAnotherGene = tmp1.isCloserToAnotherProteinCodingGene;
        saveName = 'Fig6.png';
    else
        isCloserToAnotherGene = tmp1.isCloserToAnotherGene;
        saveName = 'ExtDataFig3.png';
    end
    
    fig = createMaximisedFigure(2, [0 0 20 15]);
    axes('Position', [.22, .22, .68, .60]);
    hold on;
    xLimVal = [0-.7, nTissues+.7];


    plot(xLimVal, log10(250)*[1,1], '--', 'Color', 0.5*[1,1,1]);
    text(xLimVal(2)+.1, log10(250), '250 bp', 'Color', 0.5*[1,1,1], 'FontSize', fontSize);
    plot(xLimVal, log10(20e3)*[1,1], '--', 'Color', 0.5*[1,1,1]);
    text(xLimVal(2)+.1, log10(20e3), '20 kbp', 'Color', 0.5*[1,1,1], 'FontSize', fontSize);

    
    isOK = tmp1.iTissue>1 & ~isCloserToAnotherGene;
    swarmchart(0*tmp1.iTissue(isOK), log10(tmp1.distance_thisGene(isOK)), 'MarkerEdgeColor', sColours.closeMutation, 'MarkerFaceColor', (1+sColours.closeMutation)/2);
    isOK = tmp1.iTissue>1 & isCloserToAnotherGene;
    swarmchart(0*tmp1.iTissue(isOK), log10(tmp1.distance_thisGene(isOK)), 'MarkerEdgeColor', sColours.distantMutation, 'MarkerFaceColor', (1+sColours.distantMutation)/2);

    isOK = tmp1.iTissue>1 & ~isCloserToAnotherGene;
    swarmchart(tmp1.iTissue(isOK), log10(tmp1.distance_thisGene(isOK)), 'MarkerEdgeColor', sColours.closeMutation, 'MarkerFaceColor', (1+sColours.closeMutation)/2);
    isOK = tmp1.iTissue>1 & isCloserToAnotherGene;
    swarmchart(tmp1.iTissue(isOK), log10(tmp1.distance_thisGene(isOK)), 'MarkerEdgeColor', sColours.distantMutation, 'MarkerFaceColor', (1+sColours.distantMutation)/2);

    isOK = tmp1.iTissue==1 & ~isCloserToAnotherGene;
    h2 = swarmchart(tmp1.iTissue(isOK), log10(tmp1.distance_thisGene(isOK)), 'MarkerEdgeColor', sColours.closeMutation, 'MarkerFaceColor', (1+sColours.closeMutation)/2);
    isOK = tmp1.iTissue==1 & isCloserToAnotherGene;
    h3 = swarmchart(tmp1.iTissue(isOK), log10(tmp1.distance_thisGene(isOK)), 'MarkerEdgeColor', sColours.distantMutation, 'MarkerFaceColor', (1+sColours.distantMutation)/2);

    maxVal = 6.5; yGap = maxVal/30; ylim([0,maxVal]);
    yVal1 = 1.5*yGap + maxVal;
    yVal2 = 3.0*yGap + maxVal;
    yVal3 = 4.5*yGap + maxVal;
    yVal4 = 6.0*yGap + maxVal;
    colourBasic = .5*[1,1,1];

    for iTissue = 0:nTissues
        if (iTissue == 0)
            isOK = tmp1.iTissue>1;
        else
            isOK = tmp1.iTissue==iTissue;
        end
        ViolinGit(log10(tmp1.distance_thisGene(isOK)), iTissue, 'ViolinColor', sColours.distanceBackground, 'ViolinAlpha', 0.5, 'ShowData', false);

        xVal = iTissue + .1;
        text(xVal, yVal1, sprintf('%.1f', median(tmp1.distance_thisGene(isOK))/1e3), 'HorizontalAlignment', 'right', 'FontSize', fontSize-2, 'Color', colourBasic);
        %text(xVal, yVal1, sprintf('%.0f [%.0f-%.f]', median(tmp1.distance_thisGene(isOK))/1e3, quantile(tmp1.distance_thisGene(isOK), .25)/1e3, quantile(tmp1.distance_thisGene(isOK), .75)/1e3), 'HorizontalAlignment', 'right', 'FontSize', fontSize-2, 'Color', colourBasic);
        text(xVal, yVal2, sprintf('%.0f%%', 100*mean(tmp1.distance_thisGene(isOK)<=250)), 'HorizontalAlignment', 'right', 'FontSize', fontSize-2, 'Color', colourBasic);
        text(xVal, yVal3, sprintf('%.0f%%', 100*mean(tmp1.distance_thisGene(isOK)>20e3)), 'HorizontalAlignment', 'right', 'FontSize', fontSize-2, 'Color', colourBasic);
        text(xVal, yVal4, sprintf('%.0f%%', 100*mean(isCloserToAnotherGene(isOK))), 'HorizontalAlignment', 'right', 'FontSize', fontSize-2, 'Color', sColours.distantMutation);
    end
    text(-1, yVal1, 'median distance (kbp)', 'HorizontalAlignment', 'right', 'FontSize', fontSize-2, 'Color', colourBasic);
    text(-1, yVal2, '\leq 250 bp', 'HorizontalAlignment', 'right', 'FontSize', fontSize-2, 'Color', colourBasic);
    text(-1, yVal3, '> 20 kbp', 'HorizontalAlignment', 'right', 'FontSize', fontSize-2, 'Color', colourBasic);
    text(-1, yVal4, 'differnt closest', 'HorizontalAlignment', 'right', 'FontSize', fontSize-2, 'Color', sColours.distantMutation);


    % text(xValues+.3, yVal+0*xValues, num2str(tableTissuesWithPancancer.nSamplesWGSandRNA, '%d'), 'HorizontalAlignment', 'right', 'FontSize', fontSize-2, 'Color', .5*[1,1,1]); % {'n = '} insetad of 'n = ' will keep the space in there!
    % text(0, yVal, 'WGS+RNA', 'HorizontalAlignment', 'right', 'FontSize', fontSize-2, 'Color', .5*[1,1,1]); % {'n = '} insetad of 'n = ' will keep the space in there!
    %
    % yVal = 3*yGap + maxVal;
    % text(xValues+.3, yVal+0*xValues, num2str(tableTissuesWithPancancer.nSamplesWGS, '%d'), 'HorizontalAlignment', 'right', 'FontSize', fontSize-2, 'Color', .5*[1,1,1]); % {'n = '} insetad of 'n = ' will keep the space in there!
    % text(0, yVal, 'WGS', 'HorizontalAlignment', 'right', 'FontSize', fontSize-2, 'Color', .5*[1,1,1]); % {'n = '} insetad of 'n = ' will keep the space in there!




    yTickVal = 0:1:7; yTickLabelVal = {'1 bp', '10 bp', '100 bp', '1 kbp', '10 kbp', '100 kbp', '1 Mbp', '10 Mbp'}; set(gca, 'YTick', yTickVal, 'yTickLabel', yTickLabelVal);
    set(gca, 'XTick', 0:nTissues, 'XTickLabel', [{'Pan-cancer Solid'}; tableTissues_data1.tissuePrint], 'FontSize', fontSize)
    xlim(xLimVal); box off;



    %     hL = legend([h2, h3], {'This gene closest', 'Different gene closest'}, 'Location', 'SouthEast', 'FontSize', fontSize, 'EdgeColor', .8*[1,1,1]); %legend boxoff;
    %     hL.Position(1) = hL.Position(1) - hL.Position(3)/4;
    %     hL.Position(2) = hL.Position(2) + hL.Position(4)/2;
    %     % hL.Position(3) = hL.Position(3)*2;
    %     % hL.Position(4) = hL.Position(4)*1.2;

    hL = legend([h2, h3], {'{\itG} is the closest gene to {\itM}', 'A different gene is the closest gene to {\itM}'}, 'Location', 'South', 'FontSize', fontSize, 'EdgeColor', .8*[1,1,1]); %legend boxoff;
    %hL.Position(1) = hL.Position(1) - hL.Position(3)/4;
    hL.Position(2) = hL.Position(2) - hL.Position(4)*4;


    ylabel('Distance between gene {\itG} and mutation {\itM}'); % {'Distance between gene {\itG} and mutation {\itM}', 'of non-coding regulatory driver candidates'} 'Distance from high-CADD SNV to TSS'
    mySaveAs(fig, imagesPath, saveName)
end
%%
% old = false;
% if (old)
%     isOK = tableGencodeGenes.isCandidate;
%     tableGencodeGenesCandidates = tableGencodeGenes(isOK,:);
%     tableGencodeGenesCandidates = sortrows(tableGencodeGenesCandidates,{'tissuePrint', 'FDR'},'ascend');
%     %%
%     % tmp1 = tableMutations_candidate(tableMutations_candidate.isOK,:);
%     tmp2 = tableMutations_candidate(tableMutations_candidate.isOK,{'gene', 'candidateGenes', 'closestTSS_geneSymbol', 'closestTSS_distance', 'distance_thisGene', 'isCloserToAnotherGene', 'isHighCADD'});
%     %%
%     tmp1 = tableMutations_candidate;
%     mean(tmp1.isCloserToAnotherGene(tmp1.iTissue>1 & tmp1.isHighCADD & ~tmp1.isExcluded))
%     mean(tmp1.isCloserToAnotherGene(tmp1.iTissue==1 & tmp1.isHighCADD & ~tmp1.isExcluded))
%     %%
%     fig = createMaximisedFigure(1);
%     nR = 3; nC = 3;
%     for iTissue = 1:nTissues
%         subplot(nR, nC, iTissue);
%         isOK = tmp1.iTissue==iTissue & tmp1.isHighCADD & ~tmp1.isExcluded;
%         histogram(tmp1.distance_thisGene(isOK), 1:100:10000)
%         title(sprintf('%s: %.1f%%', tableTissues_data1.tissuePrint{iTissue}, 100*mean(tmp1.isCloserToAnotherGene(isOK))));
%     end
%     %%
%     isOK = tmp1.iTissue>1 & tmp1.isHighCADD & ~tmp1.isExcluded & tmp1.isCloserToAnotherGene;
%     median(tmp1.distance_thisGene(isOK))
%     isOK = tmp1.iTissue>1 & tmp1.isHighCADD & ~tmp1.isExcluded & ~tmp1.isCloserToAnotherGene;
%     median(tmp1.distance_thisGene(isOK))
%     isOK = tmp1.iTissue>1 & tmp1.isHighCADD & ~tmp1.isExcluded;
%     median(tmp1.distance_thisGene(isOK))
%     % isOK = tmp1.iTissue>1 & tmp1.isHighCADD & ~tmp1.isExcluded;
%     % median(tmp1.distance_thisGene(isOK))
%     %%
%     fig = createMaximisedFigure(2, [0 0 20 15]); hold on; xLimVal = 0.5+[0, nTissues]; fontSize = 12;
%     plot(xLimVal, log10(250)*[1,1], '--', 'Color', 0.5*[1,1,1]);
%     text(xLimVal(2)+.1, log10(250), '250 bp', 'Color', 0.5*[1,1,1], 'FontSize', fontSize);
%     plot(xLimVal, log10(20e3)*[1,1], '--', 'Color', 0.5*[1,1,1]);
%     text(xLimVal(2)+.1, log10(20e3), '20 kbp', 'Color', 0.5*[1,1,1], 'FontSize', fontSize);
%
%
%
%     isOK_basic = tmp1.isHighCADD & ~tmp1.isExcluded;
%     isOK = tmp1.iTissue>1 & ~tmp1.isCloserToAnotherGene;
%     swarmchart(tmp1.iTissue(isOK), log10(tmp1.distance_thisGene(isOK)), 'MarkerEdgeColor', sColours.closeMutation, 'MarkerFaceColor', (1+sColours.closeMutation)/2);
%     isOK = tmp1.iTissue>1 & tmp1.isCloserToAnotherGene;
%     swarmchart(tmp1.iTissue(isOK), log10(tmp1.distance_thisGene(isOK)), 'MarkerEdgeColor', sColours.distantMutation, 'MarkerFaceColor', (1+sColours.distantMutation)/2);
%     isOK = tmp1.iTissue==1 & ~tmp1.isCloserToAnotherGene;
%     h2 = swarmchart(tmp1.iTissue(isOK), log10(tmp1.distance_thisGene(isOK)), 'MarkerEdgeColor', sColours.closeMutation, 'MarkerFaceColor', (1+sColours.closeMutation)/2);
%     isOK = tmp1.iTissue==1 & tmp1.isCloserToAnotherGene;
%     h3 = swarmchart(tmp1.iTissue(isOK), log10(tmp1.distance_thisGene(isOK)), 'MarkerEdgeColor', sColours.distantMutation, 'MarkerFaceColor', (1+sColours.distantMutation)/2);
%
%     isOK = isOK_basic;
%     % h1 = boxchart(tmp1.iTissue(isOK), log10(tmp1.distance_thisGene(isOK)), 'BoxFaceColor', sColours.distanceBackground, 'BoxFaceAlpha', 0.5, 'LineWidth', 1.5);
%     for iTissue = 1:nTissues
%         ViolinGit(log10(tmp1.distance_thisGene(isOK & tmp1.iTissue==iTissue)), iTissue, 'ViolinColor', sColours.distanceBackground, 'ViolinAlpha', 0.5, 'ShowData', false);
%     end
%
%     % hL = legend([h1, h2, h3], {'Distribution of all pairs', 'This gene closest', 'Different gene closest'}, 'Location', 'SouthEast', 'FontSize', fontSize-2); legend boxoff; box off;
%     hL = legend([h2, h3], {'This gene closest', 'Different gene closest'}, 'Location', 'SouthEast', 'FontSize', fontSize); legend boxoff; box off;
%     hL.Position(1) = hL.Position(1) - hL.Position(3)/4;
%     hL.Position(2) = hL.Position(2) + hL.Position(4);
%
%     yTickVal = 0:1:7; yTickLabelVal = {'1 bp', '10 bp', '100 bp', '1 kbp', '10 kbp', '100 kbp', '1 Mbp', '10 Mbp'}; set(gca, 'YTick', yTickVal, 'yTickLabel', yTickLabelVal);
%     set(gca, 'XTick', 1:nTissues, 'XTickLabel', tableTissues_data1.tissuePrint, 'FontSize', fontSize)
%     xlim(xLimVal)
%     ylabel('Distance from high-CADD SNV to TSS');
%     mySaveAs(fig, imagesPath, 'Fig6b_subplotDraft.png')
%
%     % mean(tmp1.isCloserToAnotherGene(isOK_basic))
%     % mean(tmp1.isCloserToAnotherGene(isOK_basic & tmp1.iTissue>1))
%     %%
%     tmp2 = grpstats(tableMutations_candidate(tableMutations_candidate.iTissue>1 & tableMutations_candidate.isHighCADD & ~tableMutations_candidate.isExcluded, {'chrNumeric', 'pos1', 'alt'}), {'chrNumeric', 'pos1', 'alt'});
%     max(tmp2.GroupCount)
%     tmp2 = grpstats(tableMutations_candidate(tableMutations_candidate.iTissue>1 & tableMutations_candidate.isHighCADD & ~tableMutations_candidate.isExcluded, {'chrNumeric', 'pos1'}), {'chrNumeric', 'pos1'});
%     max(tmp2.GroupCount)
%     sum(tmp2.GroupCount>1)
%     tmp2 = grpstats(tableMutations_candidate(tableMutations_candidate.iTissue>1 & ~tableMutations_candidate.isExcluded, {'chrNumeric', 'pos1', 'alt'}), {'chrNumeric', 'pos1', 'alt'});
%     max(tmp2.GroupCount)
%     tmp2 = grpstats(tableMutations_candidate(tableMutations_candidate.iTissue>1, {'chrNumeric', 'pos1', 'alt'}), {'chrNumeric', 'pos1', 'alt'});
%     max(tmp2.GroupCount)
%     %%
%     isOK = tableGencodeGenesCandidates.isCandidateSolid;
%     tmp = tableGencodeGenesCandidates.nMutationsHighCADD.*tableGencodeGenesCandidates.pMutationsHighCADD_promoter;
%     n1 = sum(tmp(isOK))/sum(tableGencodeGenesCandidates.nMutationsHighCADD(isOK));
%     tmp = tableGencodeGenesCandidates.nMutationsHighCADD.*tableGencodeGenesCandidates.pMutationsHighCADD_distant;
%     n2 = sum(tmp(isOK))/sum(tableGencodeGenesCandidates.nMutationsHighCADD(isOK));
%     tmp = tableGencodeGenesCandidates.nMutationsHighCADD.*tableGencodeGenesCandidates.pMutationsHighCADD_isCloserToAnotherGene;
%     n3 = sum(tmp(isOK))/sum(tableGencodeGenesCandidates.nMutationsHighCADD(isOK));
%     fprintf('Solid cancers candidate regulatory driver mutations: %.0f%% in promoter, %.0f%% distal | %.0f%% closer to another TSS.\n', n1, n2, n3);
%     isOK = true | tableGencodeGenesCandidates.isCandidateSolid;
%     tmp = tableGencodeGenesCandidates.nMutationsHighCADD.*tableGencodeGenesCandidates.pMutationsHighCADD_promoter;
%     n1 = sum(tmp(isOK))/sum(tableGencodeGenesCandidates.nMutationsHighCADD(isOK));
%     tmp = tableGencodeGenesCandidates.nMutationsHighCADD.*tableGencodeGenesCandidates.pMutationsHighCADD_distant;
%     n2 = sum(tmp(isOK))/sum(tableGencodeGenesCandidates.nMutationsHighCADD(isOK));
%     tmp = tableGencodeGenesCandidates.nMutationsHighCADD.*tableGencodeGenesCandidates.pMutationsHighCADD_isCloserToAnotherGene;
%     n3 = sum(tmp(isOK))/sum(tableGencodeGenesCandidates.nMutationsHighCADD(isOK));
%     fprintf('All cancers candidate regulatory driver mutations: %.0f%% in promoter, %.0f%% distal | %.0f%% closer to another TSS.\n', n1, n2, n3);
%     % Solid cancers candidate regulatory driver mutations: 12.6% in promoter, 54.0% distal | 67.0% closer to another TSS.
%     % All cancers candidate regulatory driver mutations: 2.7% in promoter, 63.9% distal | 70.1% closer to another TSS.
%     %%
%     mean(tableGencodeGenesCandidates.pMutationsHighCADD_isCloserToAnotherGene==0)
%     mean(tableGencodeGenesCandidates.pMutationsHighCADD_isCloserToAnotherGene==100)
%     %%
%     mean(tableGencodeGenesCandidates.pMutationsHighCADD_promoter(tableGencodeGenesCandidates.isCandidateSolid)>0)
%     mean(tableGencodeGenesCandidates.pMutationsHighCADD_distant(tableGencodeGenesCandidates.isCandidateSolid)>0)
%     mean(tableGencodeGenesCandidates.pMutationsHighCADD_isCloserToAnotherGene(tableGencodeGenesCandidates.isCandidateSolid)>0)
%     %%
%     tableTissuesWithPancancer.nDriverUpregulatedGenes./tableTissuesWithPancancer.nCandidates
% end