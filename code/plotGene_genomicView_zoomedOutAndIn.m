function [x_mutatedUE_relative, colour] = plotGene_genomicView_zoomedOutAndIn(tissueName, biosampleABC, geneName, sColours, plotOnlyMutatedEnhancers, sProperties, exclusionType)
%%
if (~exist('exclusionType', 'var') || strcmp(exclusionType, 'excludePOLE_MSI'))
    suffix = '';
else
    suffix = ['_', exclusionType];
end
%%
% Created in saveForOneGeneVisualisation.m
load([sProperties.DIRECTORY_SAVE, '/oneGene/oneGene_', tissueName, '_', biosampleABC, '_', geneName, suffix], 'gene_pM', 'gene_qCombined', ...
    'gene_pos0', 'gene_pos1', 'gene_TSS', 'gene_strand', 'gene_nUEs', 'tableMutationsThisGene', 'tableUniqueEnhancers_oneGene', 'tableUE_annotations_hyperUE_oneGene');

% gunzip -c /share/hormozdiarilab/Codes/Regulatory_Elements/data/genes/GENCODE/gencode.v19.annotation.gtf.gz | awk '{if ($1 == "chr7" && $4>=148560000 && $5<=148640000) {print}}' > EZH2b.gencode.v19.annotation.gtf.txt

%%
tableEZH2 = readtable([sProperties.GENES_EXAMPLE_DIR, 'EZH2.gencode.v19.annotation.gtf.txt'], 'Delimiter', '\t'); % {'\t', ';', ' '}
tableEZH2.Properties.VariableNames = {'chr', 'source', 'featureType', 'pos0', 'pos1', 'score', 'strand', 'phase', 'tags'};
tableEZH2.pos0 = tableEZH2.pos0 - 1;
tableEZH2.gene_name = regexp(tableEZH2.tags, '(?<=gene_name ")[^";]*', 'once', 'match');
tableEZH2.gene_type = regexp(tableEZH2.tags, '(?<=gene_type ")[^";]*', 'once', 'match');
tableEZH2.transcript_id = regexp(tableEZH2.tags, '(?<=transcript_id ")[^";]*', 'once', 'match');
tableEZH2.transcript_name = regexp(tableEZH2.tags, '(?<=transcript_name ")[^";]*', 'once', 'match');
tableEZH2.transcript_type = regexp(tableEZH2.tags, '(?<=transcript_type ")[^";]*', 'once', 'match');
tableEZH2.isTagBasic = contains(tableEZH2.tags, 'tag "basic";');
tableEZH2.isCanonicalTranscript = contains(tableEZH2.transcript_id, 'ENST00000320356') | contains(tableEZH2.transcript_id, 'ENST00000652332') | contains(tableEZH2.transcript_id, 'ENST00000494652') | ...
    contains(tableEZH2.transcript_id, 'ENST00000365658') | contains(tableEZH2.transcript_id, 'ENST00000515903') | contains(tableEZH2.transcript_id, 'ENST00000364228') | contains(tableEZH2.transcript_id, 'ENST00000365484') | ...
    contains(tableEZH2.transcript_id, 'ENST00000516507') | contains(tableEZH2.transcript_id, 'ENST00000516501') | contains(tableEZH2.transcript_id, 'ENST00000286091'); % ENST00000652332 would be better
tableEZH2 = tableEZH2(tableEZH2.isCanonicalTranscript & ismember(tableEZH2.featureType,{'exon', 'CDS', 'UTR'}),:);
%%
rng(1);

xMin = min(tableUniqueEnhancers_oneGene.min_pos0) - 1e3;
xMax = max(tableUniqueEnhancers_oneGene.max_pos1) + 1e3; % 65e3 to see PDIA4 ENSG00000155660
%%
hold on;
lineWidth = 1.5;
minVal = 0;
maxVal = 1;
maxVal2 = maxVal/4;

cmap = lines(gene_nUEs);


if (~plotOnlyMutatedEnhancers)
    lstGenes = unique(tableEZH2.gene_name);
    [~, tableEZH2.iGene] = ismember(tableEZH2.gene_name, lstGenes);

    y_m = 0.9*maxVal;
    y_h = 0.1*maxVal;
    xStep = (xMax-xMin)/50;

    colour = 0*[1,1,1];
    for iGene = 1:length(lstGenes)
        lstRows = find(tableEZH2.iGene==iGene)';
        for iRow = lstRows
            if (~strcmp(tableEZH2.featureType{iRow}, 'CDS'))
                y_b = (y_m - y_h/2)*[1,1];
                y_t = (y_m + y_h/2)*[1,1];
            else
                y_b = (y_m - y_h)*[1,1];
                y_t = (y_m + y_h)*[1,1];
            end
            x2 = [tableEZH2.pos0(iRow), tableEZH2.pos1(iRow)];
            patch([x2, fliplr(x2)], [y_b, fliplr(y_t)], colour, 'EdgeColor', colour, 'LineWidth', lineWidth);    %set(h, 'FaceAlpha', 0.1);
        end
        tmp = sort([tableEZH2.pos0(lstRows); tableEZH2.pos1(lstRows)]);
        x2 = min(xMax, max(xMin, [tmp(1), tmp(end)]));
        if (diff(x2)>0)
            if (strcmp(tableEZH2.chr{lstRows(1)}, '+'))
                marker = '>'; xStart = x2(1);
            else
                marker = '<'; xStart = x2(2);
            end
            if (diff(x2)>xStep)
                xValues = [x2(1):xStep:x2(2), x2(2)];
                plot(xValues, y_m + 0*xValues, [marker, '-k'], 'MarkerSize', 3, 'MarkerFaceColor', 'k');
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


y_b = minVal*[1,1];
y_t = maxVal2*[1,1];

for jUE = 1:gene_nUEs
    colour = cmap(jUE, :);
    x2 = [tableUniqueEnhancers_oneGene.min_pos0(jUE), tableUniqueEnhancers_oneGene.max_pos1(jUE)];
    h = patch([x2, fliplr(x2)], [y_b, fliplr(y_t)], colour, 'EdgeColor', colour, 'LineWidth', lineWidth);    set(h, 'FaceAlpha', 0.1);

    plotCurve(gene_TSS, mean(x2), mean(y_t), maxVal, 1, colour, lineWidth);
end

isMutatedUE = tableUE_annotations_hyperUE_oneGene.foldChangeScoreM>5;
x_mutatedUE = [min(tableUniqueEnhancers_oneGene.min_pos0(isMutatedUE)), max(tableUniqueEnhancers_oneGene.max_pos1(isMutatedUE))];

colour = cmap(find(isMutatedUE, 1, 'first'), :);

if (plotOnlyMutatedEnhancers)
    isOK = tableMutationsThisGene.isHighCADD;
    if (sum(isOK)>0)
        xValues = tableMutationsThisGene.pos1(isOK);
        yValues_basic = linspace(0, maxVal2, length(xValues));
        yValues = yValues_basic(randperm(length(xValues)));
        plot(xValues, yValues, 'h', 'MarkerSize', 7, 'Color', sColours.mutated, 'MarkerFaceColor', (1+sColours.mutated)/2);
        plot((xValues*[1,1])', repmat([0, maxVal2], length(xValues), 1)', ':', 'Color', sColours.mutated);
    end

    margin = (x_mutatedUE(2) - x_mutatedUE(1))/100;
    xlim(x_mutatedUE + margin*[-1,1]);
    ylim([0, maxVal2]);
    xlabel(tableUniqueEnhancers_oneGene.chr{1})
else
    xlim([xMin, xMax]); ylim([minVal, maxVal+.3]);
    title(sprintf('{\\it%s} in %s\n\\fontsize{10}\\color[rgb]{0.5,0.5,0.5}{\\itp_M = %s, q = %s}', geneName, tissueName, getPValueAsTextTimes(gene_pM), getPValueAsTextTimes(gene_qCombined)));
end

set(gca, 'YColor', 'none');
xTickValue = get(gca, 'XTick'); xTickValue = xTickValue(1:2:end);
set(gca, 'XTick', xTickValue, 'XTickLabel', arrayfun(@num2sepNumStr, xTickValue, 'UniformOutput', false), 'XColor', 0.5*[1,1,1]);

xLimVal = get(gca, 'XLim');
x_mutatedUE_relative = (x_mutatedUE-xLimVal(1))/(xLimVal(2)-xLimVal(1));
