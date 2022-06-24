function plotSupFigure2(imagesPath, sColours, dataDepMap)


fontSize = 12;
fontSizeSmaller = fontSize - 3;
%%
fig = createMaximisedFigure(6, [0 0 30 30]); fontSize = 12;
nR = 4; nC = 2; iS = 1; xS = 0.8; yS = 0.45; xB = 0.1; yB = 0.08; xM = -0.05; yM = -0.06;

positionVector = myGeneralSubplot(nR,nC,iS,1.1+xS,yS-.05,xB,yB,xM,yM, false); iS = iS + nC; 
plotDepMap_boxplots(positionVector, 1);
positionVector = myGeneralSubplot(nR,nC,iS,1.1+xS,yS-.05,xB,yB,xM,yM, false); iS = iS + nC; 
plotDepMap_boxplots(positionVector, 2);

myGeneralSubplot(nR,nC,iS,xS,yS,xB,yB,xM,yM); iS = iS + 1; hold on;
plotDepMap_barPlot();

myGeneralSubplot(nR,nC,iS,xS,yS,xB,yB,xM,yM); iS = iS + 1; hold on;
plotDepMap_crossTissue(); 

% iS = iS + nC;
myGeneralSubplot(nR,nC,iS,.95+xS,.7,xB,yB,xM,yM); iS = iS + nC; 
plotDepMap_heatmapBreast(); 

% axPos1 = get(gca, 'Position');
% xa = (axPos1(1) + axPos1(3)*0.9333)*[1,1];
% ya = [axPos1(2), axPos1(2) + axPos1(4)];
% annotation('line',xa,ya, 'Color', 'k', 'LineWidth', 2);

% positionVector = myGeneralSubplot(nR,nC,iS,.9+xS,yS,xB,yB,xM,yM, false); iS = iS + nC; 
% plotDepMap_heatmapBoxplot(dataDepMap, positionVector); 

fontSizeLetters = 26;
dim = [.007 .99 .01 .01]; str = 'a'; annotation('textbox',dim,'String',str, 'FontSize', fontSizeLetters, 'EdgeColor','none', 'FontWeight','bold');
dim = [.007 .75 .01 .01]; str = 'b'; annotation('textbox',dim,'String',str, 'FontSize', fontSizeLetters, 'EdgeColor','none', 'FontWeight','bold');
dim = [.007 .50 .01 .01]; str = 'c'; annotation('textbox',dim,'String',str, 'FontSize', fontSizeLetters, 'EdgeColor','none', 'FontWeight','bold');
dim = [.500 .50 .01 .01]; str = 'd'; annotation('textbox',dim,'String',str, 'FontSize', fontSizeLetters, 'EdgeColor','none', 'FontWeight','bold');
dim = [.007 .27 .01 .01]; str = 'e'; annotation('textbox',dim,'String',str, 'FontSize', fontSizeLetters, 'EdgeColor','none', 'FontWeight','bold');


mySaveAs(fig, imagesPath, 'ExtDataFig2.png', false, true);
%%
    function plotDepMap_boxplots(positionVector, iTypeMetric)
        nTissues = length(dataDepMap.sResults);
        xWidth = positionVector(3)/nTissues;
        for iTissue = 1:nTissues
            labels = dataDepMap.lstGenesDepMap;
            if (iTissue == 1)
                groups = dataDepMap.groups;
                titleName = 'Pan-cancer Solid';
                if (iTypeMetric == 1)
                    yValues = dataDepMap.vDepMapGenes_average; yLabelText = 'Average gene dependency';
                else
                    yValues = dataDepMap.vDepMapGenes_pAbove50PerGene; yLabelText = 'Dependent cell lines (%)';
                end
            else
                groups = dataDepMap.sResults{iTissue}.groups;
                titleName = dataDepMap.sResults{iTissue}.tissuePrint;
                if (iTypeMetric == 1)
                    yValues = dataDepMap.matDepMatGenesTissues_average(:,iTissue);
                else
                    yValues = dataDepMap.matDepMatGenesTissues_pAbove50(:,iTissue);
                end
            end
            axes('Position', [positionVector(1)+(iTissue-1)*xWidth, positionVector(2), 0.7*xWidth, positionVector(4)]); hold on;
            
            %groups(isnan(yValues)) = NaN;
            %xValues = myJitter(groups, yValues);
            xValues = groups; 
            tmp_yValues = yValues; tmp_yValues(groups ~= 2) = NaN;
            jitter = randn(size(groups, 1), 1)/16;
            [~, perm] = sort(tmp_yValues, 'descend'); 
            lstElements = perm(1:2:end); xValues(lstElements) = xValues(lstElements) + abs(jitter(lstElements));
            lstElements = perm(2:2:end); xValues(lstElements) = xValues(lstElements) - abs(jitter(lstElements));
            

            cmap = [sColours.WT; sColours.mutated];
            for iGroup = 1:2
                if (iGroup == 2)
                    plot(xValues(groups==iGroup), yValues(groups==iGroup), 'ok', 'MarkerFaceColor', (1+cmap(iGroup,:))/2, 'Color', cmap(iGroup,:), 'MarkerSize', 4);
                end
                boxchart(iGroup*ones(sum(groups==iGroup), 1), yValues(groups==iGroup), 'BoxFaceColor', (1+cmap(iGroup,:))/2, 'LineWidth', 2, 'MarkerStyle', '+', 'MarkerColor', .6*[1,1,1], 'WhiskerLineColor', (1+cmap(iGroup,:))/2); %  % , 'Color', cmap(iGroup,:)
            end

            %boxchart(groups, yValues, 'BoxFaceColor', .5*[1,1,1], 'LineWidth', 2, 'MarkerStyle', '+', 'MarkerColor', .5*[1,1,1], 'WhiskerLineColor', .5*[1,1,1]); grid on;
            
            if (iTissue > 1)
                isOK = groups==2;
                xShift = .07;
                minValue = 0.15;
                if (iTissue == 3 && iTypeMetric == 2)
                    minValue = 7;
                end
                isOK1 = isOK & yValues > minValue & xValues > 2;
                text(xValues(isOK1) + xShift, yValues(isOK1), labels(isOK1), 'HorizontalAlignment', 'left', 'FontSize', fontSizeSmaller, 'FontAngle', 'italic');
                isOK1 = isOK & yValues > minValue & xValues <= 2;
                text(xValues(isOK1) - xShift, yValues(isOK1), labels(isOK1), 'HorizontalAlignment', 'right', 'FontSize', fontSizeSmaller, 'FontAngle', 'italic');
                set(gca, 'YTick', [])
            end
            xlim([.6,2.5]);
            set(gca, 'XTick', 1:2, 'XTickLabel', {'Control', 'Driver-up'}, 'XTickLabelRotation', 45, 'FontSize', fontSize); %, 'FontSize', 16); 
            if (iTissue == 1)
                ylabel(yLabelText);
                nCellLines = size(dataDepMap.tableCellLines, 1);
            else
                nCellLines = sum(dataDepMap.tableCellLines.iTissue==iTissue);
            end
            p = ranksum(yValues(groups==1), yValues(groups==2), 'tail', 'both');
            title(sprintf('%s\n\\fontsize{10}\\color[rgb]{0.5,0.5,0.5}{\\itp = %s}\n\\color[rgb]{0.5,0.5,0.5}{\\itn = %d}\n', titleName, getPValueAsTextTimes(p), nCellLines));
            drawnow;
        end
    end
%%
    function plotDepMap_heatmapBreast()
        iTissue = 3;
        lstBreastCancerCandidates = dataDepMap.sResults{iTissue}.geneName(dataDepMap.sResults{iTissue}.isCandidate & dataDepMap.sResults{iTissue}.isUP);
        lstRows = dataDepMap.lstGenesDepMap(dataDepMap.isCandidateGene); 
        isRowBreastHit = ismember(lstRows, lstBreastCancerCandidates);
        isColBreastCellLine = dataDepMap.tableCellLines.iTissue == iTissue & ~strcmp(dataDepMap.tableCellLines.cell_line_name, '');
        isColOther =  dataDepMap.tableCellLines.iTissue ~= iTissue;
        matCases = dataDepMap.matDepMatCandidatesCellLines(isRowBreastHit, isColBreastCellLine);
        matControl = dataDepMap.matDepMatCandidatesCellLines(isRowBreastHit, isColOther);
        nRows = size(matCases, 1);
        yLabels = lstRows(isRowBreastHit);
        xLabels = dataDepMap.tableCellLines.cell_line_name(isColBreastCellLine); 
        
        
        for iRow = 1:nRows
            try
            p = ranksum(matControl(iRow,:), matCases(iRow,:));
            pLeftTail = ranksum(matControl(iRow,:), matCases(iRow,:), 'tail', 'left');
            if (pLeftTail < 0.05)
                yLabels{iRow} = sprintf('%s %s', getPValueStarsAsText(p), yLabels{iRow});
            end
            catch
                warning('Issue in %s: %d non-missing control values, %d non-missing case values\n', yLabels{iRow}, sum(~isnan(matControl(iRow,:))), sum(~isnan(matCases(iRow,:))));
            end
        end

        [~, permCols] = sort(mean(matCases, 1, 'omitnan'), 'descend');
        [~, permRows] = sort(max(matCases, [], 2, 'omitnan'), 'descend');
        matCases = matCases(permRows, :);
        matCases = matCases(:, permCols);
        xLabels = xLabels(permCols);
        yLabels = yLabels(permRows);
        heatmap(xLabels, yLabels, matCases);
        cmap = flipud(lbmap(121, 'RedBlue')); 
        cmap = cmap([1:50, 61, 72:121],:);
        colormap(cmap); caxis([0,1]);
    end
%%
    function plotDepMap_crossTissue()
        matEnrichmentIterationTissue_crossTissue = dataDepMap.matEnrichmentIterationTissue_crossTissue(:,2:end);
        tableTissuesWithPancancer_DepMap = dataDepMap.tableTissuesWithPancancer_DepMap(2:end,:);
        nRows = size(matEnrichmentIterationTissue_crossTissue, 1);
        nTissues = size(tableTissuesWithPancancer_DepMap, 1);
        xValues = 1:nTissues;
        for iTissue = 1:nTissues
            h1 = boxchart(repmat(iTissue, nRows, 1), matEnrichmentIterationTissue_crossTissue(:,iTissue), 'BoxFaceColor', sColours.crossTissueColour, 'LineWidth', 1.5, 'MarkerStyle', '+', 'MarkerColor', sColours.crossTissueColour, 'WhiskerLineColor', sColours.crossTissueColour);
        end
        h2 = plot(xValues, tableTissuesWithPancancer_DepMap.enrichment, 'o', 'Color', sColours.mutated, 'MarkerFaceColor', (1+sColours.mutated)/2, 'MarkerSize', 7);
        %%%%%%% TODO DELETE: %%%%%%%
        %         tableTissuesWithPancancer_DepMap.crossTissue_pPermutation = NaN*ones(nTissues, 1);
        %         for iTissue = 1:nTissues
        %             p1 = mean(matEnrichmentIterationTissue_crossTissue(:,iTissue)>=tableTissuesWithPancancer_DepMap.enrichment(iTissue));
        %             p2 = mean(matEnrichmentIterationTissue_crossTissue(:,iTissue)<=tableTissuesWithPancancer_DepMap.enrichment(iTissue));
        %             p = 2*min(p1, p2);
        %             tableTissuesWithPancancer_DepMap.crossTissue_pPermutation(iTissue) = p;
        %         end
        %%%%%%% %%%%%%%%%%%% %%%%%%%
        yVal1 = 4.5;
        text(0.5, yVal1, '{\itp}-value', 'HorizontalAlignment', 'right', 'FontSize', fontSize-2, 'Color', sColours.gray);
        text(xValues, yVal1+0*xValues, arrayfun(@getPValueAsText, tableTissuesWithPancancer_DepMap.crossTissue_pPermutation, 'UniformOutput', false), 'HorizontalAlignment', 'center', 'FontSize', fontSize-2, 'Color', .5*[1,1,1]);
        
        set(gca, 'XTick', 1:nTissues, 'XTickLabel', tableTissuesWithPancancer_DepMap.tissuePrint, 'XTickLabelRotation', 45, 'FontSize', fontSize);%, 'TickLength', [0 0]);
        ylabel({'Cancer-essential', 'gene enrichment'});
        hL = legend([h2, h1], {'Matched tissue', 'Unmatched tissue'}, 'Location', 'NorthEast', 'FontSize', fontSize-2); legend boxoff
        hL.Position(2) = hL.Position(2) + hL.Position(4)*1.8;
        box off; xlim([0, nTissues]+.5);
    end
%%
    function plotDepMap_barPlot()
        tableTissuesWithPancancer_DepMap = dataDepMap.tableTissuesWithPancancer_DepMap(2:end,:);
        matValues = [tableTissuesWithPancancer_DepMap.pGenesEssential_nonDriver, tableTissuesWithPancancer_DepMap.pGenesEssential_driverUP];
        xValues = (1:size(matValues, 1))';
        yValues = max(matValues, [], 2);
        hB = bar(matValues, 'EdgeColor', 'flat', 'FaceColor', 'flat');
        hB(1).CData = sColours.gray;
        hB(2).CData = sColours.darkRed;
        text(xValues, 10 + yValues, strcat(num2str(tableTissuesWithPancancer_DepMap.enrichment, '%.1fx')), 'HorizontalAlignment', 'center', 'FontSize', fontSize-2, 'Color', .5*[1,1,1]); % {'n = '} insetad of 'n = ' will keep the space in there!
        text(xValues, 20 + yValues, arrayfun(@getPValueStarsAsText, tableTissuesWithPancancer_DepMap.pFisher, 'UniformOutput', false), 'HorizontalAlignment', 'center', 'FontSize', fontSize-2, 'Color', .5*[1,1,1]); % {'n = '} insetad of 'n = ' will keep the space in there!
        
        maxVal = max(20 + yValues); yGap = maxVal/10;
        yVal = 0.2*yGap + maxVal;        ylim([0, yVal]);
        yVal1 = 1.5*yGap + maxVal;
        yVal2 = 3*yGap + maxVal;
        yVal3 = 4.5*yGap + maxVal;

        text(xValues+.3, yVal1+0*xValues, strcat(num2str(tableTissuesWithPancancer_DepMap.nGenesEssential_control/1e3, '%.0f'), '/', num2str(tableTissuesWithPancancer_DepMap.n_control/1e3, '%.0f')), 'HorizontalAlignment', 'right', 'FontSize', fontSize-2, 'Color',  sColours.gray); % {'n = '} insetad of 'n = ' will keep the space in there!
        text(xValues+.3, yVal2+0*xValues, strcat(num2str(tableTissuesWithPancancer_DepMap.nGenesEssential_driverUP, '%d'), '/', num2str(tableTissuesWithPancancer_DepMap.n_driverUP, '%d')), 'HorizontalAlignment', 'right', 'FontSize', fontSize-2, 'Color',  sColours.darkRed); % {'n = '} insetad of 'n = ' will keep the space in there!
        

        text(0, yVal1, 'Control \times 10^3', 'HorizontalAlignment', 'right', 'FontSize', fontSize-2, 'Color', sColours.gray);
        text(0, yVal2, 'Driver-up.', 'HorizontalAlignment', 'right', 'FontSize', fontSize-2, 'Color',  sColours.darkRed);
        text(0, yVal3, 'Genes', 'HorizontalAlignment', 'right', 'FontSize', fontSize-2, 'Color', 'k');


        set(gca, 'XTick', xValues, 'XTickLabel', strrep(tableTissuesWithPancancer_DepMap.tissuePrint, 'wo Blood', 'Solid'), 'XTickLabelRotation', 45, 'FontSize', fontSize, 'TickLength', [0 0]);
        ylabel({'Tissue-matched', 'cancer-essential', 'genes (%)'});
        hL = legend({'Control', 'Driver-up'}, 'Location', 'NorthEast', 'FontSize', fontSize-2); legend boxoff
        hL.Position(2) = hL.Position(2) + hL.Position(4)/4;
        box off; xlim([0, xValues(end)+1]);
    end
end