function [hLeg, textLeg] = plotTissueScatter(sColours, sResults, iTissue, fontSize, showMoreLabels, doPlotXLabel, doPlotYLabel, doPlotLegend, lstGenesToPrint)

if (~exist('lstGenesToPrint', 'var'))
    lstGenesToPrint = {};
end
    tissuePrint = sResults{iTissue}.tissuePrint;
    pM = sResults{iTissue}.pM;
    pE = sResults{iTissue}.pE;
    qCombined = sResults{iTissue}.qCombined;
    isUP = sResults{iTissue}.isUP;
    isCandidate = sResults{iTissue}.isCandidate;
    isDriver = sResults{iTissue}.isDriver;
    geneName = sResults{iTissue}.geneName;

    
    pM(isnan(pM)) = 1; % Not mutated genes have NaN p-value
    pE(isnan(pE)) = 1; % Not mutated genes have NaN p-value
    maxValue = 16;
    xValues = -log10(pM); xValues(xValues>maxValue) = maxValue;
    yValues = -log10(pE); yValues(yValues>maxValue) = maxValue;

    if (sum(isCandidate)==0)
        maxShown_yValue = Inf;
        maxShown_xValue = Inf;
    else
        maxShown_yValue = ceil(max([5; yValues(isCandidate)]));
        maxShown_xValue = ceil(max([4; xValues(isCandidate)]));
    end

    xValues(xValues>maxShown_xValue) = maxShown_xValue; % We show all genes, but anything above the limits will be shown at the limit
    yValues(yValues>maxShown_yValue) = maxShown_yValue;

    isWithinLimits = (yValues <= maxShown_yValue) & (xValues <= maxShown_xValue); % & isM;
    plot(xValues(isWithinLimits), yValues(isWithinLimits), '.');

    xLimVal = [0, maxShown_xValue]; 
    yLimVal = [0, maxShown_yValue]; 
    alphas = [0.001, 0.01, 0.05];
    for iAlpha = 3 %1:3
        alpha = alphas(iAlpha);
        logAlpha = -log10(alpha);
        plot(xLimVal, logAlpha*[1,1], ':k'); % , 'LineWidth', lineWidth(iAlpha) text(xLimVal(2)*1.02, logAlpha, sprintf('{\\itp=%g}', alpha), 'HorizontalAlignment', 'left', 'Color', 'k', 'FontSize', fontSize);
        plot(logAlpha*[1,1], yLimVal, ':k'); % , 'LineWidth', lineWidth(iAlpha) text(logAlpha, yLimVal(2)*1.02, sprintf('{\\itp=%g}', alpha), 'HorizontalAlignment', 'center', 'Color', 'k', 'FontSize', fontSize);
    end

    textLeg = {'up', 'down', 'CDG', 'other'}; 
    
    hLeg = zeros(length(textLeg), 1);

    hLeg(1) = plot(-1,-1, '^', 'MarkerFaceColor', sColours.candidate, 'MarkerSize', 6, 'Color', sColours.candidateEdge);
    hLeg(2) = plot(-1,-1, 'v', 'MarkerFaceColor', sColours.candidate, 'MarkerSize', 6, 'Color', sColours.candidateEdge);
    hLeg(3) = plot(-1,-1, 'o', 'MarkerFaceColor', sColours.cadidateCDG, 'MarkerSize', 6, 'Color', sColours.cadidateCDGEdge);
    hLeg(4) = plot(-1,-1, 'o', 'MarkerFaceColor', sColours.other, 'MarkerSize', 6, 'Color', sColours.otherEdge);

    isOK = isDriver & isCandidate & isUP;
    if (sum(isOK)>0), plot(xValues(isOK), yValues(isOK), '^', 'MarkerFaceColor', sColours.cadidateCDG, 'MarkerSize', 6, 'Color', sColours.cadidateCDGEdge); end
    
    isOK = isDriver & isCandidate & ~isUP;
    if (sum(isOK)>0), plot(xValues(isOK), yValues(isOK), 'v', 'MarkerFaceColor', sColours.cadidateCDG, 'MarkerSize', 6, 'Color', sColours.cadidateCDGEdge); end

    isOK = isWithinLimits & ~isDriver & isCandidate & isUP;
    if (sum(isOK)>0), plot(xValues(isOK), yValues(isOK), '^', 'MarkerFaceColor', sColours.candidate, 'MarkerSize', 6, 'Color', sColours.candidateEdge); end
    
    isOK = isWithinLimits & ~isDriver & isCandidate & ~isUP;
    if (sum(isOK)>0), plot(xValues(isOK), yValues(isOK), 'v', 'MarkerFaceColor', sColours.candidate, 'MarkerSize', 6, 'Color', sColours.candidateEdge); end

    isOK = isDriver & ~isCandidate;
    if (sum(isOK)>0), plot(xValues(isOK), yValues(isOK), 'o', 'MarkerFaceColor', sColours.otherCDG, 'MarkerSize', 4, 'Color', sColours.otherCDGEdge); end

    isOK = isWithinLimits & ~isDriver & ~isCandidate;
    if (sum(isOK)>0), plot(xValues(isOK), yValues(isOK), 'o', 'MarkerFaceColor', sColours.other, 'MarkerSize', 4, 'Color', sColours.otherEdge); end
    

    xLimVal(2) = xLimVal(2)*1.01;
    yLimVal(2) = yLimVal(2)*1.01;


    if (iTissue == 1)
        xLimVal(2) = 17;
    end

    xlim(xLimVal); ylim(yLimVal);

    %     if (iTissue == 1)
    %         isShownText = isCandidate & ((isDriver & qCombined<1e-5) | qCombined<1e-10);
    %     elseif (iTissue == 3)
    %         isShownText = isCandidate & (isDriver | qCombined<0.001);
    %     else
    %         isShownText = isCandidate & (isDriver | qCombined<0.05);
    %     end
    
    if (showMoreLabels && iTissue == 1)
        isShownText = isCandidate & ((isDriver & qCombined<1e-5) | (pE<1e-6));
    elseif (~showMoreLabels && iTissue == 1)
        isShownText = isCandidate & ((isDriver & pE<1e-6) | (pM<1e-8 & pE<1e-8)); % qCombined<1e-5 & 
    elseif ismember(iTissue, [3,6])
        isShownText = isCandidate & (isDriver | qCombined<=quantile(qCombined(isCandidate & ~isDriver),1/5));
    elseif (iTissue == 7)
        isShownText = isCandidate;
    else
        isShownText = isCandidate & (isDriver | qCombined<=quantile(qCombined(isCandidate & ~isDriver),1/2));
        %     else
        %         isShownText = isCandidate;
    end
    isShownText = isShownText | (isCandidate & ismember(geneName, lstGenesToPrint));

    labels = geneName(isShownText);

    if (showMoreLabels)
        [xValuesText, yValuesText] = labelRepelSimple(xValues(isShownText), yValues(isShownText), labels, fontSize, 0.03, 0.05, 0.015, 0.5); % 0.1, 0.05, 0.015, 0.5
    else
        [xValuesText, yValuesText] = labelRepelSimple(xValues(isShownText), yValues(isShownText), labels, fontSize, 0.03, 0.05, 0.015, 0.5); % 0.1, 0.05, 0.015, 0.5
    end

    isOK2 = isDriver(isShownText);
    isToTheRight = xValuesText>xValues(isShownText);
    isOK3 = isOK2 & isToTheRight;
    text(xValuesText(isOK3), yValuesText(isOK3), labels(isOK3), 'Color', 'r', 'FontSize', fontSize, 'HorizontalAlignment', 'left');
    isOK3 = isOK2 & ~isToTheRight;
    text(xValuesText(isOK3), yValuesText(isOK3), labels(isOK3), 'Color', 'r', 'FontSize', fontSize, 'HorizontalAlignment', 'right');
    isOK3 = ~isOK2 & isToTheRight;
    text(xValuesText(isOK3), yValuesText(isOK3), labels(isOK3), 'Color', 'k', 'FontSize', fontSize, 'HorizontalAlignment', 'left');
    isOK3 = ~isOK2 & ~isToTheRight;
    text(xValuesText(isOK3), yValuesText(isOK3), labels(isOK3), 'Color', 'k', 'FontSize', fontSize, 'HorizontalAlignment', 'right');



    xTickValue = get(gca, 'XTick');
    yTickValue = get(gca, 'YTick');

    if (xLimVal(2)>10)
        isOK_x = round(xTickValue/4)==xTickValue/4;
    elseif (xLimVal(2)>6)
        isOK_x = round(xTickValue/2)==xTickValue/2;
    else
        isOK_x = round(xTickValue)==xTickValue;
    end
    if (yLimVal(2)>10)
        isOK_y = round(yTickValue/4)==yTickValue/4;
    elseif (yLimVal(2)>6)
        isOK_y = round(yTickValue/2)==yTickValue/2;
    else
        isOK_y = round(yTickValue)==yTickValue;
    end
    

    set(gca, 'FontSize', fontSize, 'XTick', xTickValue(isOK_x), 'YTick', yTickValue(isOK_y));


    if (doPlotYLabel)
        ylabel('{\itscore_E}');
    end
    if (doPlotXLabel)
        xlabel('{\itscore_M}');
    end
    if (doPlotLegend)
        h = legend(hLeg(hLeg>0), textLeg(hLeg>0), 'Location', 'NorthEast', 'FontSize', fontSize-2);
        h.Position(2) = h.Position(2) - h.Position(4);
    end
    title(tissuePrint, 'FontSize', fontSize+4);