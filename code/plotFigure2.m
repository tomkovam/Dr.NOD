function plotFigure2(imagesPath, sColours, tableTissuesWithPancancer, tableTissues_data3, tableABC, sResCrossTissues, sResCrossCADD, sResPanCancerCrossCADD, tableMutations_candidate, lstMinCADD_PHRED, sProperties)

tableTissuesWithPancancer = tableTissuesWithPancancer([end,end-1,1:end-2],:);

nRows = size(tableTissuesWithPancancer, 1);
%%
fig = createMaximisedFigure(2, [0 0 30 30]); 
fontSize = 12;
nR = 3; nC = 2; iS = 1; xS = 0.85; yS = 0.8; xB = 0.1; yB = 0.05; xM = -0.03; yM = -0.01;

myGeneralSubplot(nR,nC,1,xS,yS,xB,yB,xM,yM); hold on;
plotObservedExpected(tableTissuesWithPancancer, sColours, nRows);
axPos1 = get(gca, 'Position');

myGeneralSubplot(nR,nC,2,xS,yS,xB,yB,xM,yM); hold on;
plotEffectOfExpression(tableTissuesWithPancancer, sColours, nRows);
axPos2 = get(gca, 'Position');

axes('Position', [axPos1(1), 0.26, axPos1(3), 0.28]);
plotCrossTissue(tableTissues_data3, tableABC, sResCrossTissues);

matCrossCADD = [sResCrossCADD.enrichmentCDG; sResPanCancerCrossCADD.enrichmentCDG(2:3,:)];
lstLabels = [tableTissues_data3.tissuePrint; {'Pan-cancer Solid'; 'Pan-cancer'}];

axes('Position', [axPos2(1), 0.26, axPos2(3), 0.28]); hold on;
plotCrossCADD(matCrossCADD, lstMinCADD_PHRED, lstLabels);

% axPos2(1)+axPos2(3) % 0.9603
x = .83; y = .2;
axes('Position', [x, 0.26+y, axPos2(1)+axPos2(3)-x, 0.28-y]); hold on;
plotCrossCADD_pancancer(sResPanCancerCrossCADD, lstMinCADD_PHRED);

axes('Position', [axPos1(1), 0.05, axPos2(1)+axPos2(3)-axPos1(1), 0.15]);
plotMutationalSignatures_dotplot(tableMutations_candidate, tableTissues_data3, sProperties);


fontSizeLetters = 26;
dim = [.007 .99 .01 .01]; str = 'a'; annotation('textbox',dim,'String',str, 'FontSize', fontSizeLetters, 'EdgeColor','none', 'FontWeight','bold');
dim = [.505 .99 .01 .01]; str = 'b'; annotation('textbox',dim,'String',str, 'FontSize', fontSizeLetters, 'EdgeColor','none', 'FontWeight','bold');
dim = [.007 .58 .01 .01]; str = 'c'; annotation('textbox',dim,'String',str, 'FontSize', fontSizeLetters, 'EdgeColor','none', 'FontWeight','bold');
dim = [.505 .58 .01 .01]; str = 'd'; annotation('textbox',dim,'String',str, 'FontSize', fontSizeLetters, 'EdgeColor','none', 'FontWeight','bold');
dim = [.007 .23 .01 .01]; str = 'e'; annotation('textbox',dim,'String',str, 'FontSize', fontSizeLetters, 'EdgeColor','none', 'FontWeight','bold');



mySaveAs(fig, imagesPath, 'Fig2', false, true);
%%
    function plotCrossCADD_pancancer(sResPanCancerCrossCADD, lstMinCADD_PHRED)
        iRow = 3; typeName = 'Pan-cancer'; 
        %iRow = 2; typeName = 'Pan-cancer Solid';
        nBinCADD = find(lstMinCADD_PHRED==18);
        vecToPlot = (sResPanCancerCrossCADD.enrichmentCDG(iRow,1:nBinCADD))';
        yValues = log2(vecToPlot);
        yValues(isinf(yValues)) = NaN;

        xValues = lstMinCADD_PHRED(1:nBinCADD);
        xLimVal = [lstMinCADD_PHRED(1)-.5,lstMinCADD_PHRED(nBinCADD)+.5];
        plot(xLimVal, log2([1,1]), '-k', 'LineWidth', 2);
        
        colour = .1*[1,1,1];
        plot(xValues, yValues, 'o', 'MarkerSize', 5, 'MarkerFaceColor', (1+colour)/2, 'Color', colour);

        linearModel = fitlm(xValues, yValues); 
        x = [lstMinCADD_PHRED(1); lstMinCADD_PHRED(nBinCADD)];
        y = predict(linearModel,x);
        plot(x, y, 'k', 'LineWidth', 2, 'Color', colour);
        
        [r,p] = corr(xValues', yValues, 'type', 'Pearson');
        xLimVal = get(gca, 'XLim');
        yLimVal = get(gca, 'YLim');
        text(xLimVal(2) - .05*diff(xLimVal), yLimVal(1) + 0.1*diff(yLimVal), sprintf('{\\itr = %.1f}\n{\\itp = %s}', r, getPValueAsTextTimes(p)), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'Color', 0.5*[1,1,1])
        %         if (sum(~isnan(yValues))>2)
        %             p = linearModel.coefTest;
        %         else
        %             p = NaN;
        %         end
        fprintf('%s: r = %.1g, p = %s\n', typeName, r, getPValueAsTextShort(p));
        %fprintf('{\\bf%s}: {\\itslope=%.1f, p=%s}\n', 'Pan-cancer Solid', linearModel.Coefficients.Estimate(2), getPValueAsTextShort(p));
        title(typeName);
    end
%%
    function plotCrossCADD(matCrossCADD, lstMinCADD_PHRED, lstLabels)
        nTissues = length(lstLabels);
        %lstTissues = 1:nTissues;
        hLeg = zeros(nTissues, 1); %legendValues = cell(nTissues, 1);
        cmapTissues = linspecer(nTissues); cmapTissues(end-1,:) = cmapTissues(end,:); cmapTissues(end,:) = .1*[1,1,1];
        nBinCADD = find(lstMinCADD_PHRED==18); %length(lstMinCADD_PHRED)-2;
        %         xTestName = 'M_fullModel';
        %         yTestName = 'E_poisson';
        %         matToPlot = (sResCrossCADD.(xTestName).(yTestName).enrichmentCDG(:,1:nBinCADD))';
        matToPlot_log = log2(matCrossCADD(:,1:nBinCADD)');
        matToPlot_log(isinf(matToPlot_log)) = NaN;

        xValues = lstMinCADD_PHRED(1:nBinCADD);
        xLimVal = [lstMinCADD_PHRED(1)-.5,lstMinCADD_PHRED(nBinCADD)+.5];
        plot(xLimVal, log2([1,1]), '-k', 'LineWidth', 2);

        for iTissue = 1:size(matToPlot_log, 2)
            yValues = matToPlot_log(:,iTissue);
            colour = cmapTissues((iTissue),:);
            hLeg(iTissue) = plot(xValues, yValues, 'o', 'MarkerSize', 5, 'MarkerFaceColor', (1+colour)/2, 'Color', colour);

            %try
                linearModel = fitlm(xValues, yValues);
                x = [lstMinCADD_PHRED(1); lstMinCADD_PHRED(nBinCADD)];
                y = predict(linearModel,x);
                plot(x, y, 'k', 'LineWidth', 2, 'Color', colour);

                [r,p] = corr(xValues', yValues, 'type', 'Pearson');
                fprintf('%s: r = %.1g, p = %s\n', lstLabels{iTissue}, r, getPValueAsTextShort(p));
        %                 if (sum(~isnan(yValues))>2)
        %                     p = linearModel.coefTest;
        %                 else
        %                     p = NaN;
        %                 end
        %                 fprintf('{\\bf%s}: {\\itslope=%.1f, p=%s}\n', lstLabels{iTissue}, linearModel.Coefficients.Estimate(2), getPValueAsTextShort(p));
                drawnow
        end
        legend(hLeg, lstLabels, 'Location', 'SouthEastOutside');
        xlim(xLimVal);  grid on; grid minor;
        set(gca, 'XTick', xValues, 'XTickLabel', xValues, 'FontSize', fontSize); ylabel('Log_2 CDG enrichment');  xlabel('CADD PHRED cut-off');
    end
% CADD: 0:2:18
% {\bfBlood}: {\itslope=0.1, p=6e-10 (***)}
% {\bfBrain}: {\itslope=0.2, p=7e-06 (***)}
% {\bfBreast}: {\itslope=0.1, p=0.03 (*)}
% {\bfColorectal}: {\itslope=0.1, p=4e-05 (***)}
% {\bfLiver}: {\itslope=0.0, p=0.04 (*)}
% {\bfLung}: {\itslope=0.0, p=0.2}
% {\bfPancreas}: {\itslope=0.2, p=0.09}
% {\bfOvary}: {\itslope=0.1, p=0.07}
% CADD: 0:1:20
% {\bfBlood}: {\itslope=0.1, p=7e-05 (***)}
% {\bfBrain}: {\itslope=0.2, p=0.002 (**)}
% {\bfBreast}: {\itslope=0.1, p=0.2}
% {\bfColorectal}: {\itslope=0.1, p=0.03 (*)}
% {\bfLiver}: {\itslope=0.1, p=0.03 (*)}
% {\bfLung}: {\itslope=0.1, p=0.3}
% {\bfPancreas}: {\itslope=0.1, p=0.3}
% {\bfOvary}: {\itslope=0.1, p=0.2}
% CADD: 0:1:18
% {\bfBlood}: {\itslope=0.1, p=6e-08 (***)}
% {\bfBrain}: {\itslope=0.2, p=7e-06 (***)}
% {\bfBreast}: {\itslope=0.1, p=0.03 (*)}
% {\bfColorectal}: {\itslope=0.1, p=4e-04 (***)}
% {\bfLiver}: {\itslope=0.0, p=0.04 (*)}
% {\bfLung}: {\itslope=0.0, p=0.2}
% {\bfPancreas}: {\itslope=0.2, p=0.09}
% {\bfOvary}: {\itslope=0.1, p=0.06}
%%
    function plotCrossTissue(tableTissues, tableABC, sResTissues)
        nTissues = size(tableTissues, 1);
        %nABC = size(tableABC, 1);

        matPValue = sResTissues.pFisherCDG;
        matPValue(matPValue(:)>1) = 1;
        matPValue_log10 = -log10(matPValue);
        %matNormalized = round(matPValue_log10./max(matPValue_log10, [], 2), 2);
        %matNormalized = matPValue_log10./max(matPValue_log10, [], 2);
        matNormalized = matPValue_log10 - min(matPValue_log10, [], 2);
        matNormalized = matNormalized./max(matNormalized, [], 2);
        matToPlot = matNormalized;
        
        % Raw values normlized within row
        barTitleText = 'Row-normalized p-value';
        cmap = flipud(lbmap(300, 'RedBlue')); 
        %cmap = myColourGradient([1,1,1], sColours.darkRed, 100);
        colormap(cmap);

        imagesc(matToPlot); h = colorbar;
        set(h, 'Ticks', [0, 1], 'FontSize', fontSize);
        h.Label.String = barTitleText;
        h.Label.FontSize = fontSize;

        for iTissue = 1:nTissues
            lstABC = find(tableABC.iTissue == iTissue);
            text(mean(lstABC), 0.18, tableTissues.tissuePrint{iTissue}, 'HorizontalAlignment','left', 'FontSize', fontSize, 'Rotation', 45);
            text(0.2, iTissue, tableTissues.tissuePrint{iTissue}, 'HorizontalAlignment','right', 'FontSize', fontSize);
        end
        xlabel('Enhancer Data'); ylabel('Cancer Data'); set(gca, 'YAxisLocation', 'right')
        set(gca, 'XTick', []);
        set(gca, 'YTick', []);
        set(gca, 'FontSize', fontSize, 'XTickLabelRotation', 45, 'TickLength', [0,0]);
    end
%%
    function plotObservedExpected(tableTissuesWithPancancer, sColours, nRows)
        matValues = [tableTissuesWithPancancer.nExpected_candidateDrivers_CDGs, tableTissuesWithPancancer.nObserved_candidateDrivers_CDGs];
        xValues = (1:size(tableTissuesWithPancancer, 1))';
        yValues = max(matValues, [], 2);
        hB = bar(matValues, 'EdgeColor', 'flat', 'FaceColor', 'flat');
        hB(1).CData = sColours.gray;
        hB(2).CData = sColours.darkRed; 
        %text(xValues, 3 + yValues, tableTissuesWithPancancer.pFisherCDG_text, 'HorizontalAlignment', 'center', 'FontSize', 14);
        text(xValues, 3 + yValues, strcat(num2str(round(tableTissuesWithPancancer.enrichmentCDG), '%dx')), 'HorizontalAlignment', 'center', 'FontSize', fontSize-2, 'Color', .5*[1,1,1]); % {'n = '} insetad of 'n = ' will keep the space in there!
        text(xValues, 5 + yValues, arrayfun(@getPValueStarsAsText, tableTissuesWithPancancer.pFisherCDG, 'UniformOutput', false), 'HorizontalAlignment', 'center', 'FontSize', fontSize-2, 'Color', .5*[1,1,1]); % {'n = '} insetad of 'n = ' will keep the space in there!
        % text(xValues, 2 + yValues, strcat({'FC: '}, num2str(tableTissuesWithPancancer.enrichmentCDG, '%-.1f')), 'HorizontalAlignment', 'center', 'FontSize', 14, 'Color', .5*[1,1,1]); % {'n = '} insetad of 'n = ' will keep the space in there!
        %text(xValues, 1 + yValues, strcat({'n = '}, num2str(tableTissuesWithPancancer.nSamplesWGSandRNA, '%-d')), 'HorizontalAlignment', 'center', 'FontSize', fontSize-8, 'Color', .5*[1,1,1]); % {'n = '} insetad of 'n = ' will keep the space in there!

        maxVal = max(5 + yValues); yGap = maxVal/20;
        yVal = 0.2*yGap + maxVal;        ylim([0, yVal]);
        yVal = 1.5*yGap + maxVal;
        text(xValues+.3, yVal+0*xValues, num2str(tableTissuesWithPancancer.nSamplesWGSandRNA, '%d'), 'HorizontalAlignment', 'right', 'FontSize', fontSize-2, 'Color', .5*[1,1,1]); % {'n = '} insetad of 'n = ' will keep the space in there!
        text(0, yVal, 'WGS+RNA', 'HorizontalAlignment', 'right', 'FontSize', fontSize-2, 'Color', .5*[1,1,1]); % {'n = '} insetad of 'n = ' will keep the space in there!

        yVal = 3*yGap + maxVal;
        text(xValues+.3, yVal+0*xValues, num2str(tableTissuesWithPancancer.nSamplesWGS, '%d'), 'HorizontalAlignment', 'right', 'FontSize', fontSize-2, 'Color', .5*[1,1,1]); % {'n = '} insetad of 'n = ' will keep the space in there!
        text(0, yVal, 'WGS', 'HorizontalAlignment', 'right', 'FontSize', fontSize-2, 'Color', .5*[1,1,1]); % {'n = '} insetad of 'n = ' will keep the space in there!


        set(gca, 'XTick', xValues, 'XTickLabel', strrep(tableTissuesWithPancancer.tissuePrint, 'wo Blood', 'Solid'), 'XTickLabelRotation', 45, 'FontSize', fontSize, 'TickLength', [0 0]);
        ylabel('CDGs within targets'); %({'Cancer driver genes in targets', 'of regulatory driver candidates'});
        legend({'Expected', 'Observed'}, 'Location', 'NorthEast', 'FontSize', fontSize-2); legend boxoff %title(strrep(yTestName, '_', ' '));
        box off; xlim([0, nRows+1]);
    end
%%
    function plotEffectOfExpression(tableTissuesWithPancancer, sColours, nRows)

        tmp_tableTissues = tableTissuesWithPancancer;
        iType = 1;
        if (iType == 1)
            matValues = [tmp_tableTissues.enrichmentCDG, tmp_tableTissues.onlyP_M_enrichmentCDG];
            starsValues1 = arrayfun(@getPValueStarsAsText, tmp_tableTissues.pFisherCDG, 'UniformOutput', false);
            starsValues2 = arrayfun(@getPValueStarsAsText, tmp_tableTissues.onlyP_M_pFisherCDG, 'UniformOutput', false); 
            legendValues = {'{\itscore_M} & {\itscore_E}', '{\itscore_M} only'};
        else
            matValues = [tmp_tableTissues.pE_FDR_M_enrichmentCDG, tmp_tableTissues.onlyFDR_M_enrichmentCDG];
            legendValues = {'{\itq_M} < 2.5 & {\itp_E} < 0.05', '{\itq_M} < 2.5'};
        end
        xValues = (1:size(tmp_tableTissues, 1))';
        yValues = max(matValues, [], 2);

        plot([0, nRows+1], [1,1], ':', 'Color', .5*[1,1,1], 'LineWidth', 2);

        b = bar(matValues, 'EdgeColor', 'flat', 'FaceColor', 'flat');
        b(1).CData = sColours.darkRed;
        b(2).CData = sColours.lightRed;
        xOffset = .15; yOffset = .3; isOK = sum(matValues, 2)>0;
        text(xValues(isOK) - xOffset, yOffset + matValues(isOK,1), num2str(matValues(isOK,1), '%.1fx'), 'HorizontalAlignment', 'center', 'FontSize', fontSize-6, 'Color', .5*[1,1,1]); % {'n = '} insetad of 'n = ' will keep the space in there!
        text(xValues(isOK) + xOffset/2, yOffset + matValues(isOK,2), num2str(matValues(isOK,2), '%.1fx'), 'HorizontalAlignment', 'left', 'FontSize', fontSize-6, 'Color', .5*[1,1,1]); % {'n = '} insetad of 'n = ' will keep the space in there!
        
        text(xValues(isOK) - xOffset, 2*yOffset + matValues(isOK,1), starsValues1, 'HorizontalAlignment', 'center', 'FontSize', fontSize-6, 'Color', .5*[1,1,1]); % {'n = '} insetad of 'n = ' will keep the space in there!
        text(xValues(isOK) + xOffset/2, 2*yOffset + matValues(isOK,2), starsValues2, 'HorizontalAlignment', 'left', 'FontSize', fontSize-6, 'Color', .5*[1,1,1]); % {'n = '} insetad of 'n = ' will keep the space in there!

        yVal = 1.1*max(yValues);
        ylim([0, yVal]);
        yVal = 1.2*max(yValues);
        text(xValues(isOK), yVal+0*xValues(isOK), num2str(matValues(isOK,1)./matValues(isOK,2), '%.1f'), 'HorizontalAlignment', 'center', 'FontSize', fontSize-2, 'Color', .5*[1,1,1]); % {'n = '} insetad of 'n = ' will keep the space in there!
        text(0, yVal, 'ratio', 'HorizontalAlignment', 'center', 'FontSize', fontSize-2, 'Color', .5*[1,1,1]); % {'n = '} insetad of 'n = ' will keep the space in there!

        ylabel({'CDG enrichment'});
        set(gca, 'XTick', xValues, 'XTickLabel', strrep(tmp_tableTissues.tissuePrint, 'wo Blood', 'Solid'), 'XTickLabelRotation', 45, 'FontSize', fontSize, 'TickLength', [0 0]);
        legend(b, legendValues, 'Location', 'NorthWest', 'FontSize', fontSize-4); legend boxoff %title(strrep(yTestName, '_', ' '));
        box off; xlim([0, nRows+1]); drawnow;
    end
end