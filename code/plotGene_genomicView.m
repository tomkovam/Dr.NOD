function plotGene_genomicView(tissueName, biosampleABC, geneName, sColours, plotOnlyMutatedEnhancers, printTissue, sProperties, exclusionType)


if (~exist('exclusionType', 'var') || strcmp(exclusionType, 'excludePOLE_MSI'))
    suffix = '';
else
    suffix = ['_', exclusionType];
end
%%
% Created in saveForOneGeneVisualisation.m
load([sProperties.DIRECTORY_SAVE, '/oneGene/oneGene_', tissueName, '_', biosampleABC, '_', geneName, suffix], 'gene_pM', 'expressionPerSample', 'sampleGroup',  ...
    'gene_pos0', 'gene_pos1', 'gene_TSS', 'gene_strand', 'gene_nUEs', 'tableMutationsThisGene', 'isMutPerEnhancer', 'tableUniqueEnhancers_oneGene', 'tableUE_annotations_hyperUE_oneGene');

%%
% gunzip -c /share/hormozdiarilab/Codes/Regulatory_Elements/data/genes/GENCODE/gencode.v19.annotation.gtf.gz | awk '{if ($1 == "chr8" && $5>=128600000 && $4<=129400000) {print}}' > data/genes/MYC.gencode.v19.annotation.gtf.txt
% gunzip -c /share/hormozdiarilab/Codes/Regulatory_Elements/data/genes/GENCODE/gencode.v19.annotation.gtf.gz | awk '{if ($1 == "chr6" && $5>=132700000 && $4<=134500000) {print}}' > data/genes/SGK1.gencode.v19.annotation.gtf.txt
% gunzip -c /share/hormozdiarilab/Codes/Regulatory_Elements/data/genes/GENCODE/gencode.v19.annotation.gtf.gz | awk '{if ($1 == "chr3" && $5>=16300000 && $4<=16590000) {print}}' > data/genes/RFTN1.gencode.v19.annotation.gtf.txt
% gunzip -c /share/hormozdiarilab/Codes/Regulatory_Elements/data/genes/GENCODE/gencode.v19.annotation.gtf.gz | awk '{if ($1 == "chr6" && $5>=36940000 && $4<=37141000) {print}}' > data/genes/PIM1.gencode.v19.annotation.gtf.txt
% gunzip -c /share/hormozdiarilab/Codes/Regulatory_Elements/data/genes/GENCODE/gencode.v19.annotation.gtf.gz | grep 'tag "basic";' > data/genes/tag.basic.gencode.v19.annotation.gtf.txt

inFileOneGene = [sProperties.GENES_EXAMPLE_DIR,geneName,'.gencode.v19.annotation.gtf.extended.txt'];
if (exist(inFileOneGene, 'file'))
    inFile = inFileOneGene;
else
    inFile = [sProperties.GENES_EXAMPLE_DIR,'tag.basic.gencode.v19.annotation.gtf.txt'];
end
tableGenes = readtable(inFile, 'Delimiter', '\t'); % {'\t', ';', ' '}
lstCols = {'chr', 'source', 'featureType', 'pos0', 'pos1', 'score', 'strand', 'phase', 'tags'};
tableGenes.Properties.VariableNames = lstCols;
tableGenes.pos0 = tableGenes.pos0 - 1;
tableGenes.gene_name = regexp(tableGenes.tags, '(?<=gene_name ")[^";]*', 'once', 'match');
tableGenes.gene_type = regexp(tableGenes.tags, '(?<=gene_type ")[^";]*', 'once', 'match');
tableGenes.transcript_id = regexp(tableGenes.tags, '(?<=transcript_id ")[^";]*', 'once', 'match');
tableGenes.transcript_name = regexp(tableGenes.tags, '(?<=transcript_name ")[^";]*', 'once', 'match');
tableGenes.transcript_type = regexp(tableGenes.tags, '(?<=transcript_type ")[^";]*', 'once', 'match');
tableGenes.isTagBasic = contains(tableGenes.tags, 'tag "basic";');
tableGenes.isCanonicalTranscript = tableGenes.isTagBasic; 
tableGenes = tableGenes(tableGenes.isCanonicalTranscript & ismember(tableGenes.featureType,{'CDS', 'UTR'}),:);
minPos = min(tableUniqueEnhancers_oneGene.min_pos0);
maxPos = max(tableUniqueEnhancers_oneGene.max_pos1);
chr = tableUniqueEnhancers_oneGene.chr{1};
isOK = strcmp(tableGenes.gene_name, geneName) | (strcmp(tableGenes.chr, chr) & tableGenes.pos0 <= maxPos & tableGenes.pos1 >= minPos);
tableGenes = tableGenes(isOK,:);
if (~exist(inFileOneGene, 'file'))
    writetable(tableGenes(:,lstCols), inFileOneGene, 'Delimiter', '\t');
end
%%
fontSize = 12;
fontSizeSmaller = fontSize - 4;
hold on;

tableMutationsThisGene.yValues = tableMutationsThisGene.expression;
yValues = expressionPerSample;
maxVal_enhancer = max(yValues)*1.1;
maxVal = max(yValues)*1.5;
minVal = 0; %min(yValues)*0.9;
cmap = lines(gene_nUEs);


lstMarkers = {'o', 's', 'd', 'h'};
%%
% tableUE_annotations_hyperUE_oneGene(:,{'name', 'foldChangeScoreM'})
% tableUniqueEnhancers_oneGene
if (plotOnlyMutatedEnhancers)
    minFC = 2^10;
    if (strcmp(geneName, 'RFTN1'))
        minFC = 2^7;
    end
    if (strcmp(geneName, 'BCL2'))
        minFC = 1e7;
    end
    if (strcmp(geneName, 'EBF1'))
        minFC = 50;
    end
    if (ismember(geneName, {'HIST1H2BG', ''}))
        minFC = 0;
    end
    isMutatedUE = tableUE_annotations_hyperUE_oneGene.foldChangeScoreM>minFC;
    %tableUE_annotations_hyperUE_oneGene(isMutatedUE,:)
    x2 = [min(tableUniqueEnhancers_oneGene.min_pos0(isMutatedUE)), max(tableUniqueEnhancers_oneGene.max_pos1(isMutatedUE))];
    margin = (x2(2) - x2(1))/50;
    xLimVal = x2 + margin*[-1,1];
    %
    lstGenes = unique(tableGenes.gene_name);
    [~, tableGenes.iGene] = ismember(tableGenes.gene_name, lstGenes);
    
    
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



hLeg = []; hLegText = {'TSS'};

y_b = minVal*[1,1];
y_t = maxVal_enhancer*[1,1];

for jUE = 1:gene_nUEs
    colour = cmap(jUE, :);
    isOK = isMutPerEnhancer(:,jUE);
    plot(tableMutationsThisGene.pos1(isOK), tableMutationsThisGene.yValues(isOK), lstMarkers{2}, 'MarkerFaceColor', colour, 'Color', (colour+1)/2);

    x2 = [tableUniqueEnhancers_oneGene.min_pos0(jUE), tableUniqueEnhancers_oneGene.max_pos1(jUE)];
    h = patch([x2, fliplr(x2)], [y_b, fliplr(y_t)], colour, 'EdgeColor', colour);    set(h, 'FaceAlpha', 0.1);
    
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
    jUE = find(isMutPerEnhancer(iMut,:));
    if (length(jUE)>1)
        error('One mutation in multiple enhancers!');
    end
    h = plot(tableMutationsThisGene.pos1(iMut)*[1,1], yLimVal, ':', 'Color', cmap(jUE,:));
end
if (exist('h', 'var'))
    hLeg = [hLeg, h]; hLegText = [hLegText, {'sample without expression'}];
end

ylim(yLimVal);

if (plotOnlyMutatedEnhancers)
    xlim(xLimVal);
end
xTickValue = get(gca, 'XTick'); xTickValue = xTickValue(1:2:end);
set(gca, 'XTick', xTickValue, 'XTickLabel', arrayfun(@num2sepNumStr, xTickValue, 'UniformOutput', false), 'XColor', 0.5*[1,1,1]);



h = plot(get(gca, 'XLim'), median(yValues(sampleGroup == 1), 'omitnan')*[1,1], '--', 'Color', sColours.WT);
hLeg = [hLeg, h]; hLegText = [hLegText, {'median WT expression'}];


ax.XAxis.FontSize = fontSizeSmaller;
ax.YAxis.FontSize = fontSize;
ylabel('Expression'); xlabel(sprintf('Chr%d', tableMutationsThisGene.chrNumeric(1))); % , 'FontSize', fontSize

if (printTissue)
    titleText = sprintf('%s in %s', geneName, tissueName);
else
    titleText = geneName;
end
title(sprintf('%s\n\\fontsize{8}\\color[rgb]{0.5,0.5,0.5}{\\itp_M = %s}', titleText, getPValueAsTextTimes(gene_pM)), 'FontSize', fontSize);

