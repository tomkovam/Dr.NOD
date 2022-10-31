function plotFigure6(tableMutationGenePairs, imagesPath, sColours, tableTissues_data1)

rng(1);
tmp1 = tableMutationGenePairs(tableMutationGenePairs.isHighCADD & ~tableMutationGenePairs.isExcluded, :);

nTissues = size(tableTissues_data1, 1);
fontSize = 12;

for iType = 1:2
    if (iType == 1)
        isCloserToAnotherGene = tmp1.isCloserToAnotherProteinCodingGene;
        saveName = 'Fig7';
    else
        isCloserToAnotherGene = tmp1.isCloserToAnotherGene;
        saveName = 'ExtDataFig3';
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
    mySwarmchart(0*tmp1.iTissue(isOK), log10(tmp1.distance_thisGene(isOK)), sColours.closeMutation, (1+sColours.closeMutation)/2);
    isOK = tmp1.iTissue>1 & isCloserToAnotherGene;
    mySwarmchart(0*tmp1.iTissue(isOK), log10(tmp1.distance_thisGene(isOK)), sColours.distantMutation, (1+sColours.distantMutation)/2);

    isOK = tmp1.iTissue>1 & ~isCloserToAnotherGene;
    mySwarmchart(tmp1.iTissue(isOK), log10(tmp1.distance_thisGene(isOK)), sColours.closeMutation, (1+sColours.closeMutation)/2);
    isOK = tmp1.iTissue>1 & isCloserToAnotherGene;
    mySwarmchart(tmp1.iTissue(isOK), log10(tmp1.distance_thisGene(isOK)), sColours.distantMutation, (1+sColours.distantMutation)/2);

    isOK = tmp1.iTissue==1 & ~isCloserToAnotherGene;
    h2 = mySwarmchart(tmp1.iTissue(isOK), log10(tmp1.distance_thisGene(isOK)), sColours.closeMutation, (1+sColours.closeMutation)/2);
    isOK = tmp1.iTissue==1 & isCloserToAnotherGene;
    h3 = mySwarmchart(tmp1.iTissue(isOK), log10(tmp1.distance_thisGene(isOK)), sColours.distantMutation, (1+sColours.distantMutation)/2);

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
    text(-1, yVal4, 'different closest', 'HorizontalAlignment', 'right', 'FontSize', fontSize-2, 'Color', sColours.distantMutation);


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
    mySaveAs(fig, imagesPath, saveName, true, true);
    savefig([imagesPath, saveName, '.fig']);
end