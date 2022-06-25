function plotBarsWithLiterature(sColours, tableGencodeGenes)

yValueText = 'sizeEffectE';
tmpToPlot = tableGencodeGenes(tableGencodeGenes.isCandidateSolid,:);
tmpToPlot = sortrows(tmpToPlot,yValueText,'ascend');


xValues = (1:size(tmpToPlot, 1))';
yValues = tmpToPlot.(yValueText);

hold on; 

alphaPrognostic = 0.05; 
alphaPrognosticStrict = 0.01;
alphaPrognosticMostStrict = 0.001;

labelsGeneNames = tmpToPlot.geneSymbol;
labelsTissue = tmpToPlot.tissuePrint;

labelsTissue(abs(yValues)<0.5 & strcmp(labelsTissue, 'Colorectal')) = {'Color.'};
labelsTissue(abs(yValues)<0.5 & strcmp(labelsTissue, 'Pancreas')) = {'Panc.'};
labelsTissue(abs(yValues)<0.4 & strcmp(labelsTissue, 'Breast')) = {'Bre.'};


hB = bar(yValues);
hB.EdgeColor = 'none';
hB.FaceColor = 'flat';
hB.CData(:) = repmat(sColours.nonCDG, size(tmpToPlot, 1), 1);

hLeg = zeros(8,1); iLeg = 1;

verbose = false;
if (verbose)
    isUP = yValues>0;
    isDOWN = yValues<0;
    fprintf('%d/%d (%.1f%%) upregulated genes have oncogenic tissue-specific evidence.\n', sum(tmpToPlot.literatureEvidenceOncogene(isUP)>0), sum(isUP), 100*mean(tmpToPlot.literatureEvidenceOncogene(isUP)>0));
    fprintf('%d/%d (%.1f%%) downregulated genes have TSG tissue-specific evidence.\n', sum(tmpToPlot.literatureEvidenceTSG(isDOWN)>0), sum(isDOWN), 100*mean(tmpToPlot.literatureEvidenceTSG(isDOWN)>0));
    fprintf('%d/%d (%.1f%%) upregulated genes have strong oncogenic tissue-specific evidence.\n', sum(tmpToPlot.literatureEvidenceOncogene(isUP)>2), sum(isUP), 100*mean(tmpToPlot.literatureEvidenceOncogene(isUP)>2));
    fprintf('%d/%d (%.1f%%) downregulated genes have strong TSG tissue-specific evidence.\n', sum(tmpToPlot.literatureEvidenceTSG(isDOWN)>2), sum(isDOWN), 100*mean(tmpToPlot.literatureEvidenceTSG(isDOWN)>2));
    for tissue = unique(tmpToPlot.tissuePrint)'
        isOK = strcmp(tmpToPlot.tissuePrint, tissue{1});
        isOK2 = isOK & isUP & tmpToPlot.literatureEvidenceOncogene>2;
        fprintf('UP %s: %d genes, %d strong oncogenes (%.0f%%): %s\n', tissue{1}, sum(isOK & isUP), sum(isOK2), 100*sum(isOK2)/sum(isOK), strjoin(unique(tmpToPlot.geneSymbol(isOK2)), ', '));
    end
end


for iRow = find(tmpToPlot.literatureEvidenceOncogene>0 | tmpToPlot.literatureEvidenceTSG>0)'   
    levelONC_orig = tmpToPlot.literatureEvidenceOncogene(iRow); 
    levelTSG_orig = tmpToPlot.literatureEvidenceTSG(iRow);

    levelONC = levelONC_orig;
    levelONC(ismember(levelONC_orig, [1,2])) = 1;
    levelONC(ismember(levelONC_orig, [3,4])) = 2;

    levelTSG = levelTSG_orig;
    levelTSG(ismember(levelTSG_orig, [1,2])) = 1;
    levelTSG(ismember(levelTSG_orig, [3,4])) = 2;

    if (levelTSG==2)
        disp(levelTSG)
    end

    maxLevel = 2;    
    levelTop = NaN;
    levelBottom = NaN;

    if (levelONC>0)
        white = (maxLevel-levelONC);
        colourTop = (sColours.ONCOGENE + white)/(1+white);
        levelTop = levelONC;
    end
    if (levelTSG>0)
        white = (maxLevel-levelTSG);
        colourBottom = (sColours.TSG + white)/(1+white);
        levelBottom = maxLevel + levelTSG;
    end

    if (levelONC==0)
        colourTop = colourBottom;
        levelTop = NaN;
    end
    if (levelTSG==0)
        colourBottom = colourTop;
        levelBottom = NaN;
    end

    x2 = iRow + .4*[-1,1];
    
    if (yValues(iRow)>0)
        y_b = [0,0];
        y_t = [0, yValues(iRow)];
    else
        y_t = [yValues(iRow), 0];
        y_b = yValues(iRow)*[1,1];
    end
    h = patch([x2, fliplr(x2)], [y_b, fliplr(y_t)], colourBottom, 'EdgeColor', 'none');    %set(h, 'FaceAlpha', 0.5);
    if (~isnan(levelBottom)), hLeg(levelBottom) = h; end

    if (yValues(iRow)>0)
        y_b = [0, yValues(iRow)];
        y_t = yValues(iRow)*[1,1];
    else
        y_t = [0,0];
        y_b = [yValues(iRow), 0];
    end
    h = patch([x2, fliplr(x2)], [y_b, fliplr(y_t)], colourTop, 'EdgeColor', 'none');    %set(h, 'FaceAlpha', 0.5);
    if (~isnan(levelTop)), hLeg(levelTop) = h; end
end


legValues = {'oncogene weak evidence', 'oncogene strong evidence', 'TSG weak evidence', 'TSG strong evidence'}; 
iLeg = 5;
xVal = max(xValues) + 2;
yStep = max(abs(yValues))/30;


%%
yValuesStar = 0 - 2 * sign(yValues) * yStep; markerSizeStar = 6;
% yValuesStar = maxY + 2*stepY + 0*yValues;

isOK = tmpToPlot.isONCOGENE & ~tmpToPlot.isTSG; marker = 'o';
hLeg(iLeg) = plot(xValues(isOK), yValuesStar(isOK), marker, 'Color', sColours.ONCOGENE, 'MarkerFaceColor', 'none', 'MarkerSize', markerSizeStar, 'LineWidth', 1.5);
% legValues{iLeg} = 'CGC oncogene'; iLeg = iLeg + 1;

isOK = tmpToPlot.isTSG & ~tmpToPlot.isONCOGENE;
hLeg(iLeg) = plot(xValues(isOK), yValuesStar(isOK), marker, 'Color', sColours.TSG, 'MarkerFaceColor', 'none', 'MarkerSize', markerSizeStar, 'LineWidth', 1.5);
% legValues{iLeg} = 'CGC TSG'; iLeg = iLeg + 1;

isOK = tmpToPlot.isONCOGENE & tmpToPlot.isTSG;
plot(xValues(isOK)+.15, yValuesStar(isOK), marker, 'Color', sColours.TSG, 'MarkerFaceColor', 'none', 'MarkerSize', markerSizeStar, 'LineWidth', 1.5);
plot(xValues(isOK)-.15, yValuesStar(isOK), marker, 'Color', sColours.ONCOGENE, 'MarkerFaceColor', 'none', 'MarkerSize', markerSizeStar, 'LineWidth', 1.5);

isOK = tmpToPlot.isInCGC & ~tmpToPlot.isONCOGENE & ~tmpToPlot.isTSG;
hLeg(iLeg) = plot(xValues(isOK), yValuesStar(isOK), marker, 'Color', sColours.fusion, 'MarkerFaceColor', 'none', 'MarkerSize', markerSizeStar, 'LineWidth', 1.5);
% legValues{iLeg} = 'CGC fusion'; iLeg = iLeg + 1;

yVal = yValuesStar(end);
text(xVal, yVal, ['CGC: ', insertColourIntoText('oncogene', sColours.ONCOGENE), ' | ', insertColourIntoText('TSG', sColours.TSG), ' | ', insertColourIntoText('fusion', sColours.fusion)]);

%%
yValuesStar = 0 - 5 * sign(yValues) * yStep;
charF = char(9829); % empty heart: 9825
charU = char(10013); % greek cross: 10010); %128327); %10015); %10013); % skull: 9760

sColours.TSG_light = (1+sColours.TSG)/2;
sColours.TSG_veryLight = (4+sColours.TSG)/5;
sColours.ONCOGENE_light = (1+sColours.ONCOGENE)/2;
sColours.ONCOGENE_veryLight = (4+sColours.ONCOGENE)/5;

tmpToPlot.iClassSurvival = zeros(size(tmpToPlot, 1), 1);
tmpToPlot.iClassSurvival(tmpToPlot.pValuePrognostic_tissueMatched<alphaPrognostic & tmpToPlot.isFavorable_tissueMatched==1) = 1;
tmpToPlot.iClassSurvival(tmpToPlot.pValuePrognostic_tissueMatched<alphaPrognosticStrict & tmpToPlot.isFavorable_tissueMatched==1) = 2;
tmpToPlot.iClassSurvival(tmpToPlot.pValuePrognostic_tissueMatched<alphaPrognosticMostStrict & tmpToPlot.isFavorable_tissueMatched==1) = 3;
tmpToPlot.iClassSurvival(tmpToPlot.pValuePrognostic_tissueMatched<alphaPrognostic & tmpToPlot.isFavorable_tissueMatched==0) = 4;
tmpToPlot.iClassSurvival(tmpToPlot.pValuePrognostic_tissueMatched<alphaPrognosticStrict & tmpToPlot.isFavorable_tissueMatched==0) = 5;
tmpToPlot.iClassSurvival(tmpToPlot.pValuePrognostic_tissueMatched<alphaPrognosticMostStrict & tmpToPlot.isFavorable_tissueMatched==0) = 6;

% grpstats(tmpToPlot(:,'iClassSurvival'), 'iClassSurvival')

isOK = tmpToPlot.iClassSurvival == 1;
text(xValues(isOK), yValuesStar(isOK), charF, 'Color', sColours.TSG_veryLight, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
isOK = tmpToPlot.iClassSurvival == 2;
text(xValues(isOK), yValuesStar(isOK), charF, 'Color', sColours.TSG_light, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
isOK = tmpToPlot.iClassSurvival == 3;
text(xValues(isOK), yValuesStar(isOK), charF, 'Color', sColours.TSG, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');

yVal = yValuesStar(end);
text(xVal, yVal, [insertColourIntoText(['favourable ', charF], sColours.TSG), ': ', insertColourIntoText(sprintf('p<%.1g', alphaPrognostic), sColours.TSG_veryLight), ...
    ' | ', insertColourIntoText(sprintf('p<%.1g', alphaPrognosticStrict), sColours.TSG_light), ' | ', insertColourIntoText(sprintf('p<%.1g', alphaPrognosticMostStrict), sColours.TSG)]);

isOK = tmpToPlot.iClassSurvival == 4;
text(xValues(isOK), yValuesStar(isOK), charU, 'Color', sColours.ONCOGENE_veryLight, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
isOK = tmpToPlot.iClassSurvival == 5;
text(xValues(isOK), yValuesStar(isOK), charU, 'Color', sColours.ONCOGENE_light, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
isOK = tmpToPlot.iClassSurvival == 6;
text(xValues(isOK), yValuesStar(isOK), charU, 'Color', sColours.ONCOGENE, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');

yValuesStar = 0 - 8 * sign(yValues) * yStep;
yVal = yValuesStar(end);
text(xVal, yVal, [insertColourIntoText(['unfavourable ', charU], sColours.ONCOGENE), ': ', insertColourIntoText(sprintf('p<%.1g', alphaPrognostic), sColours.ONCOGENE_veryLight), ...
    ' | ', insertColourIntoText(sprintf('p<%.1g', alphaPrognosticStrict), sColours.ONCOGENE_light), ' | ', insertColourIntoText(sprintf('p<%.1g', alphaPrognosticMostStrict), sColours.ONCOGENE)]);

%%
fontSizeGenes = 10;
yValuesText = yValues + sign(yValues) * max(abs(yValues))/50;
isOK = yValues < 0;
text(xValues(isOK), yValuesText(isOK), labelsGeneNames(isOK), 'Rotation', 90, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', 'FontSize', fontSizeGenes, 'FontAngle', 'italic');
isOK = yValues > 0;
text(xValues(isOK), yValuesText(isOK), labelsGeneNames(isOK), 'Rotation', 90, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'FontSize', fontSizeGenes, 'FontAngle', 'italic');


fontSizeTissues = 8;
yValuesText = sign(yValues) * max(abs(yValues))/50;
isOK = yValues > 0;
text(xValues(isOK), yValuesText(isOK), labelsTissue(isOK), 'Color', 'w', 'Rotation', 90, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'FontSize', fontSizeTissues);
isOK = yValues < 0;
text(xValues(isOK), yValuesText(isOK), labelsTissue(isOK), 'Color', 'w', 'Rotation', 90, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', 'FontSize', fontSizeTissues);
%%
fontSizeLabels = 14;
set(gca, 'FontSize', fontSizeLabels, 'XTick', [], 'XColor', 'none', 'TickLength', [0 0]);

ylabel('Expression size effect');

lstOK = 1:4;
hLegend = legend(hLeg(lstOK), legValues(lstOK), 'Location', 'EastOutside', 'FontSize', 10);
title(hLegend,'Tissue-specific evidence');
hLegend.Title.Visible = 'on';
