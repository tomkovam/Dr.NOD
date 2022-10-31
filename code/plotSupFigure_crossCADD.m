function plotSupFigure_crossCADD(imagesPath, tableTissues_data3, sResCrossCADD, lstMinCADD_PHRED)


nTissues = size(sResCrossCADD.enrichmentCDG, 1);

fig = createMaximisedFigure(1, [0 0 30 15]);

nR = 2; nC = 4; iS = 1; xS = 0.75; yS = 0.75; xB = 0.05; yB = 0.1; xM = -0.02; yM = -0.02;

fontSize = 10;

for iTissue = 1:nTissues
    %subplot(nR,nC,iTissue);
    myGeneralSubplot(nR,nC,iTissue,xS,yS,xB,yB,xM,yM); hold on;

    matCrossCADD = sResCrossCADD.enrichmentCDG(iTissue,:);
    plotCrossCADD(matCrossCADD, lstMinCADD_PHRED);
    title(tableTissues_data3.tissuePrint{iTissue});
end
mySaveAs(fig, imagesPath, 'SupFig_crossCADD.png', true, true);
savefig([imagesPath, 'SupFig_crossCADD.fig']);
%%
    function plotCrossCADD(matCrossCADD, lstMinCADD_PHRED)
        nBinCADD = find(lstMinCADD_PHRED==18);
        matToPlot_log = log2(matCrossCADD(:,1:nBinCADD)');
        matToPlot_log(isinf(matToPlot_log)) = NaN;

        xValues = lstMinCADD_PHRED(1:nBinCADD);
        xLimVal = [lstMinCADD_PHRED(1)-.5,lstMinCADD_PHRED(nBinCADD)+.5];
        plot(xLimVal, log2([1,1]), '-k', 'LineWidth', 2);

        yValues = matToPlot_log;
        colour = 0*[1,1,1]; 
        plot(xValues, yValues, 'o', 'MarkerSize', 5, 'MarkerFaceColor', (1+colour)/2, 'Color', colour);

        linearModel = fitlm(xValues, yValues);
        x = [lstMinCADD_PHRED(1); lstMinCADD_PHRED(nBinCADD)];
        y = predict(linearModel,x);
        plot(x, y, 'k', 'LineWidth', 2, 'Color', colour);

        [r,p] = corr(xValues', yValues, 'type', 'Pearson', 'rows', 'pairwise');
        fprintf('%s: r = %.1g, p = %s\n', tableTissues_data3.tissuePrint{iTissue}, r, getPValueAsTextShort(p));
        drawnow

        xlim(xLimVal); ylim([-4,4]); grid on; grid minor;
        set(gca, 'XTick', xValues, 'XTickLabel', xValues, 'FontSize', fontSize, 'XTickLabelRotation', 0); 
        ylabel('Log_2 CDG enrichment');  xlabel('CADD PHRED cut-off');
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
%     function plotCrossCADD_pancancer(sResPanCancerCrossCADD, lstMinCADD_PHRED)
%         iRow = 3; typeName = 'Pan-cancer';
%         %iRow = 2; typeName = 'Pan-cancer Solid';
%         nBinCADD = find(lstMinCADD_PHRED==18);
%         vecToPlot = (sResPanCancerCrossCADD.enrichmentCDG(iRow,1:nBinCADD))';
%         yValues = log2(vecToPlot);
%         yValues(isinf(yValues)) = NaN;
% 
%         xValues = lstMinCADD_PHRED(1:nBinCADD);
%         xLimVal = [lstMinCADD_PHRED(1)-.5,lstMinCADD_PHRED(nBinCADD)+.5];
%         plot(xLimVal, log2([1,1]), '-k', 'LineWidth', 2);
% 
%         colour = .1*[1,1,1];
%         plot(xValues, yValues, 'o', 'MarkerSize', 5, 'MarkerFaceColor', (1+colour)/2, 'Color', colour);
% 
%         linearModel = fitlm(xValues, yValues);
%         x = [lstMinCADD_PHRED(1); lstMinCADD_PHRED(nBinCADD)];
%         y = predict(linearModel,x);
%         plot(x, y, 'k', 'LineWidth', 2, 'Color', colour);
% 
%         [r,p] = corr(xValues', yValues, 'type', 'Pearson');
%         xLimVal = get(gca, 'XLim');
%         yLimVal = get(gca, 'YLim');
%         text(xLimVal(2) - .05*diff(xLimVal), yLimVal(1) + 0.1*diff(yLimVal), sprintf('{\\itr = %.1f}\n{\\itp = %s}', r, getPValueAsTextTimes(p)), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'Color', 0.5*[1,1,1])
%         %         if (sum(~isnan(yValues))>2)
%         %             p = linearModel.coefTest;
%         %         else
%         %             p = NaN;
%         %         end
%         fprintf('%s: r = %.1g, p = %s\n', typeName, r, getPValueAsTextShort(p));
%         %fprintf('{\\bf%s}: {\\itslope=%.1f, p=%s}\n', 'Pan-cancer Solid', linearModel.Coefficients.Estimate(2), getPValueAsTextShort(p));
%         title(typeName);
%     end
end