function plotGene_genomicView(tissueName, biosampleABC, geneName, sColours, plotOnlyMutatedEnhancers, printTissue, sProperties)

% Created in saveForOneGeneVisualisation.m
load(['save/oneGene/oneGene_', tissueName, '_', biosampleABC, '_', geneName], 'gene_pM', 'expressionPerSample', 'sampleGroup',  ...
    'gene_pos0', 'gene_pos1', 'gene_TSS', 'gene_strand', 'gene_nUEs', 'tableMutationsThisGene', 'isMutPerEnhancer', 'tableUniqueEnhancers_oneGene', 'tableUE_annotations_hyperUE_oneGene');
% load(['save/oneGene_', tissueName, biosampleABC, geneName], 'gene_pM', 'gene_pE', 'gene_FDR', 'isMutatedSample', 'expressionPerSample', 'CNVperSample', 'nSamples', 'sampleGroup', 'sampleGroupInclWoExpression', ...
%     'gene_pos0', 'gene_pos1', 'gene_TSS', 'gene_strand', 'gene_nUEs', 'tableMutationsThisGene', 'isMutPerEnhancer', 'tableUniqueEnhancers_oneGene', 'tableUE_annotations_hyperUE_oneGene');

%%
% gunzip -c /share/hormozdiarilab/Codes/Regulatory_Elements/data/genes/GENCODE/gencode.v19.annotation.gtf.gz | awk '{if ($1 == "chr8" && $5>=128600000 && $4<=129400000) {print}}' > data/genes/MYC.gencode.v19.annotation.gtf.txt
% gunzip -c /share/hormozdiarilab/Codes/Regulatory_Elements/data/genes/GENCODE/gencode.v19.annotation.gtf.gz | awk '{if ($1 == "chr6" && $5>=132700000 && $4<=134500000) {print}}' > data/genes/SGK1.gencode.v19.annotation.gtf.txt
% gunzip -c /share/hormozdiarilab/Codes/Regulatory_Elements/data/genes/GENCODE/gencode.v19.annotation.gtf.gz | awk '{if ($1 == "chr3" && $5>=16300000 && $4<=16590000) {print}}' > data/genes/RFTN1.gencode.v19.annotation.gtf.txt
% gunzip -c /share/hormozdiarilab/Codes/Regulatory_Elements/data/genes/GENCODE/gencode.v19.annotation.gtf.gz | awk '{if ($1 == "chr6" && $5>=36940000 && $4<=37141000) {print}}' > data/genes/PIM1.gencode.v19.annotation.gtf.txt
tableGenes = readtable([sProperties.GENES_EXAMPLE_DIR,geneName,'.gencode.v19.annotation.gtf.txt'], 'Delimiter', '\t'); % {'\t', ';', ' '}
tableGenes.Properties.VariableNames = {'chr', 'source', 'featureType', 'pos0', 'pos1', 'score', 'strand', 'phase', 'tags'};
tableGenes.pos0 = tableGenes.pos0 - 1;
tableGenes.gene_name = regexp(tableGenes.tags, '(?<=gene_name ")[^";]*', 'once', 'match');
tableGenes.gene_type = regexp(tableGenes.tags, '(?<=gene_type ")[^";]*', 'once', 'match');
tableGenes.transcript_id = regexp(tableGenes.tags, '(?<=transcript_id ")[^";]*', 'once', 'match');
tableGenes.transcript_name = regexp(tableGenes.tags, '(?<=transcript_name ")[^";]*', 'once', 'match');
tableGenes.transcript_type = regexp(tableGenes.tags, '(?<=transcript_type ")[^";]*', 'once', 'match');
tableGenes.isTagBasic = contains(tableGenes.tags, 'tag "basic";');
tableGenes.isCanonicalTranscript = tableGenes.isTagBasic; %contains(tableGenes.transcript_id, 'ENST00000320356') | contains(tableGenes.transcript_id, 'ENST00000652332') | contains(tableGenes.transcript_id, 'ENST00000494652') | ...
%     contains(tableGenes.transcript_id, 'ENST00000365658') | contains(tableGenes.transcript_id, 'ENST00000515903') | contains(tableGenes.transcript_id, 'ENST00000364228') | contains(tableGenes.transcript_id, 'ENST00000365484') | ...
%     contains(tableGenes.transcript_id, 'ENST00000516507') | contains(tableGenes.transcript_id, 'ENST00000516501') | contains(tableGenes.transcript_id, 'ENST00000286091'); % ENST00000652332 would be better
% unique(tableGenes(tableGenes.isCanonicalTranscript, {'gene_name', 'transcript_name', 'transcript_type'})) % misc_RNA
tableGenes = tableGenes(tableGenes.isCanonicalTranscript & ismember(tableGenes.featureType,{'CDS', 'UTR'}),:);
%%
fontSize = 12;
fontSizeSmaller = fontSize - 4;
hold on;

% tableUniqueEnhancers_oneGene
% tableUE_annotations_hyperUE_oneGene
% yValues = log2(1+expressionPerSample);

tableMutationsThisGene.yValues = tableMutationsThisGene.expression;
yValues = expressionPerSample;
maxVal_enhancer = max(yValues)*1.1;
maxVal = max(yValues)*1.5;
minVal = 0; %min(yValues)*0.9;
cmap = lines(gene_nUEs);

x1 = sort([gene_pos0, gene_pos1]);

lstMarkers = {'o', 's', 'd', 'h'};
%%

% tableUE_annotations_hyperUE_oneGene
if (plotOnlyMutatedEnhancers)
    minFC = 2^10;
    if (strcmp(geneName, 'RFTN1'))
        minFC = 2^7;
    end
    isMutatedUE = tableUE_annotations_hyperUE_oneGene.foldChangeScoreM>minFC;
    x2 = [min(tableUniqueEnhancers_oneGene.min_pos0(isMutatedUE)), max(tableUniqueEnhancers_oneGene.max_pos1(isMutatedUE))];
    margin = (x2(2) - x2(1))/50;
    %margin= 1e3;
    xLimVal = x2 + margin*[-1,1];
    %
    lstGenes = unique(tableGenes.gene_name);
    [~, tableGenes.iGene] = ismember(tableGenes.gene_name, lstGenes);

    %     isMutatedUE = tableUE_annotations_hyperUE_oneGene.foldChangeScoreM>5;
    %     xMin = min(tableUniqueEnhancers_oneGene.min_pos0(isMutatedUE));
    %     xMax = max(tableUniqueEnhancers_oneGene.max_pos1(isMutatedUE));
    
    xMin = xLimVal(1);
    xMax = xLimVal(2);

    y_m = 0.9*maxVal;
    y_h = 0.1*maxVal;
    xStep = (xMax-xMin)/50;
    lineWidth = 1.5;

    colour = 0*[1,1,1];
    for iGene = 1:length(lstGenes)
        lstRows = find(tableGenes.iGene==iGene)';
        for iRow = lstRows
            if (strcmp(tableGenes.featureType{iRow}, 'UTR'))
                y_b = (y_m - y_h/2)*[1,1];
                y_t = (y_m + y_h/2)*[1,1];
            else
                y_b = (y_m - y_h)*[1,1];
                y_t = (y_m + y_h)*[1,1];
            end
            x2 = [tableGenes.pos0(iRow), tableGenes.pos1(iRow)];
            patch([x2, fliplr(x2)], [y_b, fliplr(y_t)], colour, 'EdgeColor', colour, 'LineWidth', lineWidth);    %set(h, 'FaceAlpha', 0.1);            
        end
        tmp = sort([tableGenes.pos0(lstRows); tableGenes.pos1(lstRows)]);
        x2 = min(xMax, max(xMin, [tmp(1), tmp(end)]));
        if (diff(x2)>0)
            if (strcmp(tableGenes.strand{lstRows(1)}, '+'))
                marker = '>'; xStart = x2(1);
            else
                marker = '<'; xStart = x2(2);
            end
            if (diff(x2)>xStep)
                xValues = [x2(1):xStep:x2(2), x2(2)];
                plot(xValues, y_m + 0*xValues, [marker, '-k'], 'MarkerSize', 5, 'MarkerFaceColor', 'k');
            end
            if (ismember(lstGenes{iGene}, {'EZH2', 'RNY5'}))
                alignment = 'center';
            else
                alignment = 'left';
            end
            text(xStart, y_t(1), lstGenes{iGene}, 'FontAngle', 'italic', 'HorizontalAlignment', alignment, 'VerticalAlignment', 'bottom');

        end
    end
end

% if (strcmp(lstGenes{iGene}, 'MYC'))
%     tableUE_annotations_hyperUE_oneGene.foldChangeScoreM
%     x2
% end
%


hLeg = []; hLegText = {'TSS'};
% hLeg(1) = plot(gene_TSS*[1,1], [minVal, maxVal], '-', 'LineWidth', 2, 'Color', 'y');
% if (strcmp(gene_strand, '+'))
%     text(gene_TSS, 0.95*maxVal, sprintf(' >> %s >>  >>  >>  >> ', geneName), 'Color', 'y', 'HorizontalAlignment', 'left');
% else
%     text(gene_TSS, 0.95*maxVal, sprintf(' <<  <<  <<  << %s << ', geneName), 'Color', 'y', 'HorizontalAlignment', 'right');
% end

y_b = minVal*[1,1];
y_t = maxVal_enhancer*[1,1];

for jUE = 1:gene_nUEs
    colour = cmap(jUE, :);
    %     iUE = lst_iUE(iColour);
    %     isOK = tableMutationsThisGene.iUniqueEnhancer == iUE;
    isOK = isMutPerEnhancer(:,jUE);
    plot(tableMutationsThisGene.pos1(isOK), tableMutationsThisGene.yValues(isOK), lstMarkers{2}, 'MarkerFaceColor', colour, 'Color', (colour+1)/2);

    x2 = [tableUniqueEnhancers_oneGene.min_pos0(jUE), tableUniqueEnhancers_oneGene.max_pos1(jUE)];
    h = patch([x2, fliplr(x2)], [y_b, fliplr(y_t)], colour, 'EdgeColor', colour);    set(h, 'FaceAlpha', 0.1);

    %         th = linspace( pi/2, -pi/2, 100);
    %         R = 1;  %or whatever radius you want
    %         x = R*cos(th) + 5;
    %         y = R*sin(th) + 4;
    %         plot(x,y);

    % BETTER BUT NOT COMPUTED ATM: label = {sprintf('%s, p=%s, mut=%d', tableUniqueEnhancers.name{iUE}, getPValueAsText(tableUniqueEnhancers.pBinomTestRight_SNVs_highCADD(iUE)), sum(isOK))};
    label = {sprintf('%s, mut=%d, FC:%.1fx', tableUniqueEnhancers_oneGene.name{jUE}, sum(isOK), tableUE_annotations_hyperUE_oneGene.foldChangeScoreM(jUE))};
    hLeg = [hLeg, h]; hLegText = [hLegText, label];
end


lstRows = find(strcmp(tableGenes.featureType, 'CDS') & tableGenes.pos1 >= xLimVal(1) & tableGenes.pos0 <= xLimVal(2))';
for iRow = lstRows
    y_b = 0*[1,1];
    y_t = (y_m - y_h)*[1,1]; % maxVal*[1,1];
    x2 = [tableGenes.pos0(iRow), tableGenes.pos1(iRow)];
    h = patch([x2, fliplr(x2)], [0*y_b, fliplr(y_t)], .7*[1,1,1], 'EdgeColor', 'none');  %set(h, 'FaceAlpha', 0.1);
end

isOK = tableMutationsThisGene.isHighCADD;
if (sum(isOK)>0)
    h = plot(tableMutationsThisGene.pos1(isOK), tableMutationsThisGene.yValues(isOK), lstMarkers{4}, 'MarkerSize', 10, 'Color', .2*[1,1,1]);
    hLeg = [hLeg, h]; hLegText = [hLegText, {'high CADD'}];
end
isOK = tableMutationsThisGene.isIndel;
if (sum(isOK)>0)
    h = plot(tableMutationsThisGene.pos1(isOK), tableMutationsThisGene.yValues(isOK), lstMarkers{3}, 'MarkerSize', 10, 'Color', .2*[1,1,1]);
    for iMut = find(isOK)'
        plot([tableMutationsThisGene.pos0(iMut), tableMutationsThisGene.pos1(iMut)], [tableMutationsThisGene.yValues(iMut), tableMutationsThisGene.yValues(iMut)], '-', 'MarkerSize', 10, 'Color', .2*[1,1,1]);
    end
    hLeg = [hLeg, h]; hLegText = [hLegText, {'indel'}];
end

yLimVal = [minVal, maxVal]; clear h
for iMut = find(isnan(tableMutationsThisGene.yValues))'
    % OLD: [~, jUE] = ismember(tableMutationsThisGene.iUniqueEnhancer(iMut), lst_iUE); 
    jUE = find(isMutPerEnhancer(iMut,:));
    if (length(jUE)>1)
        error('One mutation in multiple enhancers!');
    end
    %if (tmp2_mutations.isIndel(iMut)), marker = lstMarkers{3};
    %elseif (tmp2_mutations.isHighCADD(iMut)), marker = lstMarkers{4}; else, marker = lstMarkers{2}; end [marker, ':']
    h = plot(tableMutationsThisGene.pos1(iMut)*[1,1], yLimVal, ':', 'Color', cmap(jUE,:));
end
if (exist('h', 'var'))
    hLeg = [hLeg, h]; hLegText = [hLegText, {'sample without expression'}];
end

%OLD legend([{'TSS'}; tableUniqueEnhancers.name(lst_iUE); {'high CADD'; 'known driver'; 'median WT'}], 'Location', 'EastOutside');


% text(tableMutationsThisGene.pos1, tableMutationsThisGene.yValues, tableMutationsThisGene.patternName);
ylim(yLimVal);


% scaleLength = 20e3;
% %xStart = scaleLength*round(((minPos+maxPos)/2 - scaleLength/2)/scaleLength); yMiddle = yLimVal(2)*0.8; plot(xStart + [0,scaleLength], yMiddle*[1,1], '-k', 'LineWidth', 2);
% xStart = gene_TSS; yMiddle = yLimVal(2)*0.8; plot(xStart + [0,scaleLength], yMiddle*[1,1], '-k', 'LineWidth', 2);
% text(xStart, yMiddle*1.01, {sprintf('%d kbp', scaleLength/1e3), ''}, 'HorizontalAlignment', 'center');

% xTickVal = get(gca, 'XTick');
% set(gca, 'XTickLabel', strcat(num2sepNumStr(round(xTickVal/1e3)), 'k'));

if (plotOnlyMutatedEnhancers)
    xlim(xLimVal);
end
xTickValue = get(gca, 'XTick'); xTickValue = xTickValue(1:2:end);
set(gca, 'XTick', xTickValue, 'XTickLabel', arrayfun(@num2sepNumStr, xTickValue, 'UniformOutput', false), 'XColor', 0.5*[1,1,1]);



% xLimVal = get(gca, 'XLim');

% h = patch([x1, fliplr(x1)], [0.90*maxVal*[1,1], fliplr(0.92*maxVal*[1,1])], [1,1,0], 'EdgeColor', 'none');    set(h, 'FaceAlpha', 0.5);   drawnow;
%text(mean(x1), yLimVal(2)*.99, getPValueAsText(tableUniqueEnhancers.min_pBinomTestRight(iUE)), 'HorizontalAlignment', 'center', 'Rotation', 90);

% xlim(xLimVal);

h = plot(get(gca, 'XLim'), median(yValues(sampleGroup == 1), 'omitnan')*[1,1], '--', 'Color', sColours.WT);
hLeg = [hLeg, h]; hLegText = [hLegText, {'median WT expression'}];




% set(gca, 'FontSize', fontSize);
ax = gca;
ax.XAxis.FontSize = fontSizeSmaller;
ax.YAxis.FontSize = fontSize;
ylabel('Expression'); xlabel(sprintf('Chr%d', tableMutationsThisGene.chrNumeric(1))); % , 'FontSize', fontSize
% legend(hLeg, hLegText, 'Location', 'EastOutside');
% mySaveAs(fig, imagesPath, [tissueName, '_', geneName,biosampleABC,'.png']);
% title(sprintf('%s in %s\n\\color[rgb]{0.5,0.5,0.5}{\\itp_M = %s}', geneName, tissueName, getPValueAsTextShort(gene_pM)));

if (printTissue)
    titleText = sprintf('%s in %s', geneName, tissueName);
else
    titleText = geneName;
end
% title(sprintf('%s\n\\color[rgb]{0.5,0.5,0.5}{\\itp_M = %s}', titleText, getPValueAsTextShort(gene_pM)));
title(sprintf('%s\n\\fontsize{8}\\color[rgb]{0.5,0.5,0.5}{\\itp_M = %s}', titleText, getPValueAsTextTimes(gene_pM)), 'FontSize', fontSize);

