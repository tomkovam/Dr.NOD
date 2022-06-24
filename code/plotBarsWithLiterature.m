function plotBarsWithLiterature(sColours, tableGencodeGenes)

% lstTypes = {'sizeEffectExpression_woBlood', 'sizeEffectExpression_onlyBlood', 'sizeEffectExpression_withBlood'};
% iType = 1;
% yValueText = lstTypes{iType}; %'sizeEffectExpression_woBlood';
% tmpToPlot = tableGencodeGenes(tableGencodeGenes.isCandidate & ~isnan(tableGencodeGenes.(yValueText)),:);

yValueText = 'sizeEffectE';

tmpToPlot = tableGencodeGenes(tableGencodeGenes.isCandidateSolid,:);
% yValueText = 'sign_pValueExpression_woBlood';
% tmpToPlot.(yValueText) = sign(tmpToPlot.sizeEffectExpression_woBlood).* -log10(tmpToPlot.pValueExpression_woBlood);


% tmpToPlot.yValue = sign(tmpToPlot.(yValueText)).*-log10(tmpToPlot.candidate_FDR);
% yValueText = 'yValue';

tmpToPlot = sortrows(tmpToPlot,yValueText,'ascend');


xValues = (1:size(tmpToPlot, 1))';
yValues = tmpToPlot.(yValueText);
% if (sum(~isnan(yValues))>0)
%     fig = createMaximisedFigure(2, [0 0 35 15]);
%     axes('Position', [0.13 0.05 0.85 0.9]);
hold on; 
alphaPrognostic = 0.01; 
alphaPrognosticStrict = 0.001;

alphaPrognostic = 0.05; 
alphaPrognosticStrict = 0.01;
alphaPrognosticMostStrict = 0.001;

labelsGeneNames = tmpToPlot.geneSymbol;
labelsTissue = tmpToPlot.tissuePrint;

labelsTissue(abs(yValues)<0.5 & strcmp(labelsTissue, 'Colorectal')) = {'Color.'};
labelsTissue(abs(yValues)<0.5 & strcmp(labelsTissue, 'Pancreas')) = {'Panc.'};
labelsTissue(abs(yValues)<0.4 & strcmp(labelsTissue, 'Breast')) = {'Bre.'};

% if (iType == 1)
%     labelsTissue = strrep(labelsTissue, 'Blood|', ''); % Here, we do not show blood results
% end


%     hLeg = [];

hB = bar(yValues);
hB.EdgeColor = 'none';
hB.FaceColor = 'flat';
hB.CData(:) = repmat(sColours.nonCDG, size(tmpToPlot, 1), 1);
% hB.CData(tmpToPlot.isDriver,:) = repmat(sColours.fusion, sum(tmpToPlot.isDriver), 1);
% hB.CData(tmpToPlot.isTSG,:) = repmat(sColours.TSG, sum(tmpToPlot.isTSG), 1);
% hB.CData(tmpToPlot.isONCOGENE,:) = repmat(sColours.ONCOGENE, sum(tmpToPlot.isONCOGENE), 1);
% hB.CData(tmpToPlot.literatureEvidenceOncogene>0 & yValues>0,:) = repmat(sColours.ONCOGENE, sum(tmpToPlot.literatureEvidenceOncogene>0 & yValues>0), 1);
% hB.CData(tmpToPlot.literatureEvidenceTSG>0 & yValues<0,:) = repmat(sColours.TSG, sum(tmpToPlot.literatureEvidenceTSG>0 & yValues<0), 1);
% hB.CData(tmpToPlot.literatureEvidenceOncogene>tmpToPlot.literatureEvidenceTSG,:) = repmat(sColours.ONCOGENE, sum(tmpToPlot.literatureEvidenceOncogene>tmpToPlot.literatureEvidenceTSG), 1);
% hB.CData(tmpToPlot.literatureEvidenceOncogene<tmpToPlot.literatureEvidenceTSG,:) = repmat(sColours.TSG, sum(tmpToPlot.literatureEvidenceOncogene<tmpToPlot.literatureEvidenceTSG), 1);

% maxY = max(abs(yValues));
% stepY = maxY/7;

hLeg = zeros(8,1); iLeg = 1;

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
    %     levelONC = min([maxLevel,levelONC]);
    %     levelTSG = min([maxLevel,levelTSG]);
    
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
% legValues = {'oncogene weak', 'oncogene strong', 'TSG weak', 'TSG strong'}; 
iLeg = 5;
xVal = max(xValues) + 2;
yStep = max(abs(yValues))/30;

% yValuesStar = yValues + sign(yValues) * yStep;
% yValuesStar = -1*sign(yValues) * yStep;


% yValuesStar = yValues - sign(yValues) * yStep;
% yValuesStar = maxY + 3*stepY + 0*yValues;

% isOK = tmpToPlot.pValuePrognostic_tissueMatched<alphaPrognostic & tmpToPlot.isFavorable_tissueMatched==1;
% text(xValues(isOK), yValuesStar(isOK), 'h', 'Color', sColours.TSG, 'MarkerFaceColor', 'w', 'MarkerSize', markerSizeStar, 'LineWidth', 1.5);
% isOK = tmpToPlot.pValuePrognostic_tissueMatched<alphaPrognostic & tmpToPlot.isFavorable_tissueMatched==0;
% text(xValues(isOK), yValuesStar(isOK), 'h', 'Color', sColours.ONCOGENE, 'MarkerFaceColor', 'w', 'MarkerSize', markerSizeStar, 'LineWidth', 1.5);


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

% yLimVal = get(gca, 'YLim'); 
% yLimVal(2) = maxY;
% ylim(yLimVal);

fontSizeLabels = 14;
set(gca, 'FontSize', fontSizeLabels, 'XTick', [], 'XColor', 'none', 'TickLength', [0 0]);

% ax = gca;
% ax.Clipping = 'off';

%     xlabel('Target genes of regulatory driver candidates');
ylabel('Expression size effect');
%     yLimVal = get(gca, 'YLim');
%     text(-4, yLimVal(1), 'Down', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', fontSizeLabels); % {'Downregulated', 'in mutated'}
%     text(-4, yLimVal(2), 'Up', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', 'FontSize', fontSizeLabels); % {'Upregulated', 'in mutated'}

lstOK = 1:4;
hLegend = legend(hLeg(lstOK), legValues(lstOK), 'Location', 'EastOutside', 'FontSize', 10);
title(hLegend,'Tissue-specific evidence');
hLegend.Title.Visible = 'on';
% hLegend.Title.NodeChildren.Position = [0.5 1.3 0];


% legend({'a', 'b', 'c', 'd', 'e', 'f'})

% %     mySaveAs(fig, imagesPath, ['Fig3_', lstTypes{iType}], false, true);
% end
%%



%%
% tmp_UP =
%
%   13×26 table
%
%     chromosome       pos0          pos1       geneSymbol     geneType1    strand       geneNameGencode            geneType2         geneNameGencodeFirstPart    isUsedGene    isCandidate    isDriver    isUP     isDOWN    isONCOGENE    isTSG    sizeEffectExpression_woBlood    pValueExpression_woBlood    nTissues_candidate    tissuePrint    typePrognostic_tissueMatched    pValuePrognostic_tissueMatched    isFavorable_tissueMatched    indexFavorable_tissueMatched_strict    indexFavorable_tissueMatched_veryStrict          label
%     __________    __________    __________    ___________    _________    ______    ______________________    __________________    ________________________    __________    ___________    ________    _____    ______    __________    _____    ____________________________    ________________________    __________________    _________________    ____________________________    ______________________________    _________________________    ___________________________________    _______________________________________    __________________
%
%     {'chr1' }     2.3884e+07    2.3886e+07    {'ID3'    }    {'gene'}     {'-'}     {'ENSG00000117318.8' }    {'protein_coding'}      {'ENSG00000117318'}         true           true         true       true     false       false       true                1.0198                        0.013045                   1             {'Ovary'        }              {'U**' }                         0.006467                           0                                 1                                       NaN                      {'...ID3|U**'    }
%     {'chr10'}     7.5669e+07    7.5677e+07    {'PLAU'   }    {'gene'}     {'+'}     {'ENSG00000122861.11'}    {'protein_coding'}      {'ENSG00000122861'}         true           true         false      true     false       false       false              0.94792                       0.0060198                   1             {'Lung'         }              {'U***'}                        0.0004036                           0                                 1                                         1                      {'...PLAU|U***'  }
%     {'chr12'}     1.2526e+08    1.2537e+08    {'SCARB1' }    {'gene'}     {'-'}     {'ENSG00000073060.11'}    {'protein_coding'}      {'ENSG00000073060'}         true           true         false      true     false       false       false               1.0631                       2.474e-05                   1             {'Liver'        }              {'U***'}                        1.831e-05                           0                                 1                                         1                      {'...SCARB1|U***'}
%     {'chr14'}     9.2432e+07    9.2507e+07    {'TRIP11' }    {'gene'}     {'-'}     {'ENSG00000100815.8' }    {'protein_coding'}      {'ENSG00000100815'}         true           true         true       true     false       false       false               0.4939                        0.020277                   1             {'Liver'        }              {'U**' }                         0.003491                           0                                 1                                       NaN                      {'...TRIP11|U**' }
%     {'chr17'}     7.9992e+06    8.0224e+06    {'ALOXE3' }    {'gene'}     {'-'}     {'ENSG00000179148.5' }    {'protein_coding'}      {'ENSG00000179148'}         true           true         false      true     false       false       false              0.78266                        0.029784                   1             {'Lung'         }              {'U**' }                        0.0003556                           0                                 1                                         1                      {'...ALOXE3|U**' }
%     {'chr3' }     1.9331e+08    1.9342e+08    {'OPA1'   }    {'gene'}     {'+'}     {'ENSG00000198836.4' }    {'protein_coding'}      {'ENSG00000198836'}         true           true         false      true     false       false       false              0.32536                        0.043912                   1             {'Head and neck'}              {'U**' }                         0.004499                           0                                 1                                       NaN                      {'...OPA1|U**'   }
%     {'chr5' }     1.4887e+08    1.4893e+08    {'CSNK1A1'}    {'gene'}     {'-'}     {'ENSG00000113712.12'}    {'protein_coding'}      {'ENSG00000113712'}         true           true         true       true     false       false       false              0.74952                        0.024366                   1             {'Head and neck'}              {'U**' }                         0.006277                           0                                 1                                       NaN                      {'...CSNK1A1|U**'}
%     {'chr6' }     1.0518e+08    1.0531e+08    {'HACE1'  }    {'gene'}     {'-'}     {'ENSG00000085382.7' }    {'protein_coding'}      {'ENSG00000085382'}         true           true         false      true     false       false       false              0.55828                        0.046342                   1             {'Liver'        }              {'U**' }                         0.001219                           0                                 1                                       NaN                      {'...HACE1|U**'  }
%     {'chr7' }      1.485e+08    1.4858e+08    {'EZH2'   }    {'gene'}     {'-'}     {'ENSG00000106462.6' }    {'protein_coding'}      {'ENSG00000106462'}         true           true         true       true     false       false       false              0.61608                       0.0055936                   1             {'Liver'        }              {'U***'}                        2.105e-06                           0                                 1                                         1                      {'...EZH2|U***'  }
%     {'chr9' }       3.51e+07    3.5103e+07    {'STOML2' }    {'gene'}     {'-'}     {'ENSG00000165283.11'}    {'protein_coding'}      {'ENSG00000165283'}         true           true         false      true     false       false       false              0.50767                        0.015234                   1             {'Liver'        }              {'U***'}                        0.0004768                           0                                 1                                         1                      {'...STOML2|U***'}
%
%     {'chr19'}     4.5971e+07    4.5978e+07    {'FOSB'   }    {'gene'}     {'+'}     {'ENSG00000125740.9' }    {'protein_coding'}      {'ENSG00000125740'}         true           true         false      true     false       false       false               1.7318                        0.004805                   1             {'Liver'        }              {'F**' }                         0.001298                           1                                 2                                       NaN                      {'...FOSB|F**'   }
%     {'chr20'}     5.0956e+06    5.1073e+06    {'PCNA'   }    {'gene'}     {'-'}     {'ENSG00000132646.6' }    {'protein_coding'}      {'ENSG00000132646'}         true           true         false      true     false       false       false              0.84043                        0.033453                   1             {'Head and neck'}              {'F**' }                         0.007109                           1                                 2                                       NaN                      {'...PCNA|F**'   }
%     {'chr3' }     1.9545e+08    1.9547e+08    {'MUC20'  }    {'gene'}     {'+'}     {'ENSG00000176945.12'}    {'protein_coding'}      {'ENSG00000176945'}         true           true         true       true     false       false       false               1.4527                      0.00054733                   1             {'Head and neck'}              {'F**' }                         0.003978                           1                                 2                                       NaN                      {'...MUC20|F**'  }
%
%
% tmp_DOWN =
%
%   1×26 table
%
%     chromosome       pos0          pos1       geneSymbol    geneType1    strand       geneNameGencode            geneType2         geneNameGencodeFirstPart    isUsedGene    isCandidate    isDriver    isUP     isDOWN    isONCOGENE    isTSG    sizeEffectExpression_woBlood    pValueExpression_woBlood    nTissues_candidate    tissuePrint    typePrognostic_tissueMatched    pValuePrognostic_tissueMatched    isFavorable_tissueMatched    indexFavorable_tissueMatched_strict    indexFavorable_tissueMatched_veryStrict         label
%     __________    __________    __________    __________    _________    ______    ______________________    __________________    ________________________    __________    ___________    ________    _____    ______    __________    _____    ____________________________    ________________________    __________________    _________________    ____________________________    ______________________________    _________________________    ___________________________________    _______________________________________    ________________
%
%      {'chrX'}     7.0338e+07    7.0362e+07    {'MED12'}     {'gene'}     {'+'}     {'ENSG00000184634.11'}    {'protein_coding'}      {'ENSG00000184634'}         true           true         true       false    true        false       true               -0.70921                      0.027169                    1                 {'Brain'}                  {'F**'}                          0.006888                           1                                 2                                       NaN                      {'...MED12|F**'}
%% Liver has many unfavorable prognostic genes, but only a few favorable prognostic genes
% tmp2 = sortrows(grpstats(tablePrognosticGenes(strcmp(tablePrognosticGenes.Cancer, 'liver cancer'), {'typePrognosticShort', 'iTypePrognostic'}), 'typePrognosticShort', 'mean'))
%% Head&neck has many favorable prognostic genes, but only a few unfavorable prognostic genes
% tmp2 = sortrows(grpstats(tablePrognosticGenes(strcmp(tablePrognosticGenes.Cancer, 'head and neck cancer'), {'typePrognosticShort', 'iTypePrognostic'}), 'typePrognosticShort', 'mean'))
% tmp2 = sortrows(grpstats(tablePrognosticGenes(strcmp(tablePrognosticGenes.Cancer, 'thyroid cancer'), {'typePrognosticShort', 'iTypePrognostic'}), 'typePrognosticShort', 'mean'))
%%
% geneSymbol = 'ALOXE3';
% iRow = find(strcmp(tablePrognosticGenes.Cancer, tableTissues.tissuePrognostic{7}) & strcmp(tablePrognosticGenes.GeneName, geneSymbol));
% tablePrognosticGenes(iRow,:)


%%
% tmpToPlot = tableGencodeGenes(tableGencodeGenes.isCandidate,:);
% % isOK = tmpToPlot.isCandidate;% & tableGencodeGenes.isDriver;
% xValues = tmpToPlot.sizeEffectExpression_woBlood; %xValues(~isOK) = NaN;
% yValues = -log10(tmpToPlot.pValueExpression_woBlood); %yValues(~isOK) = NaN;
% labels = tmpToPlot.label;
% fig = createMaximisedFigure(2); hold on; hLeg = NaN*ones(4,1);
% hLeg(1) = plot(xValues, yValues, 'o');
% isOK = tmpToPlot.isDriver;
% hLeg(2) = plot(xValues(isOK), yValues(isOK), 'or');
% isOK = tmpToPlot.isDriver & ~tmpToPlot.isONCOGENE & ~tmpToPlot.isTSG;
% text(xValues(isOK), yValues(isOK), labels(isOK), 'Color', (1+[1,0,0])/2);
% % isOK = tableGencodeGenes.isCandidate;
% % text(xValues(isOK), yValues(isOK), tableGencodeGenes.geneSymbol(isOK));
% isOK = tmpToPlot.isONCOGENE;
% hLeg(3) = plot(xValues(isOK), yValues(isOK), 'or', 'MarkerFaceColor', 'r');
% text(xValues(isOK), yValues(isOK), labels(isOK), 'Color', 'r');
% isOK = tmpToPlot.isTSG;
% hLeg(4) = plot(xValues(isOK), yValues(isOK), 'og', 'MarkerFaceColor', 'g');
% text(xValues(isOK), yValues(isOK), labels(isOK), 'Color', 'g');
% isOK = ~isnan(tmpToPlot.indexFavorable_tissueMatched_veryStrict) & ~tmpToPlot.isDriver;
% text(xValues(isOK), yValues(isOK), labels(isOK), 'Color', 'b');
% isOK = yValues > 4 & isnan(tmpToPlot.indexFavorable_tissueMatched_veryStrict) & ~tmpToPlot.isDriver;
% text(xValues(isOK), yValues(isOK), labels(isOK), 'Color', .5*[1,1,1]);
% legend(hLeg, {'not CDG', 'other CDG', 'oncogene CDG', 'TSG CDG'});
% tmpToPlot(isOK,:)
% isOK =  ~isnan(tmpToPlot.indexFavorable_tissueMatched_strict) & tmpToPlot.isUP;
% tmp_UP = tmpToPlot(isOK,:);
% tmp_UP = sortrows(tmp_UP,'isFavorable_tissueMatched','ascend') % 10 correct, 4 wrong
% isOK = ~isnan(tmpToPlot.indexFavorable_tissueMatched_strict) & tmpToPlot.isDOWN;
% tmp_DOWN = tmpToPlot(isOK,:);
% tmp_DOWN = sortrows(tmp_DOWN,'isFavorable_tissueMatched','ascend') % 1 correct, 0 wrong