function [xAltRelative, yAltRelative] = plotGene_boxplot_forLogos(tissueName, biosampleABC, geneName, sColours, iSampleMutated, exclusionType)


if (~exist('exclusionType', 'var') || strcmp(exclusionType, 'excludePOLE_MSI'))
    suffix = '';
else
    suffix = ['_', exclusionType];
end

% Created in saveForOneGeneVisualisation.m
load(['save/oneGene/oneGene_', tissueName, '_', biosampleABC, '_', geneName, suffix], 'gene_pM', 'gene_pE', 'gene_qCombined', 'expressionPerSample', 'nSamples', 'sampleGroup', 'sampleGroupInclWoExpression');

fontSize = 10;

xLimVal = 1 + [-0.4, 0.5];

hold on;
xValues = 1+randn(nSamples, 1)/16;
yValues = expressionPerSample; yLabelText = 'FPKM_{uq}';

cmap = [sColours.WT; sColours.mutated];
for iGroup = 1
    plot(xValues(sampleGroup==iGroup), yValues(sampleGroup==iGroup), 'ok', 'MarkerFaceColor', (1+cmap(iGroup,:))/2, 'Color', cmap(iGroup,:), 'MarkerSize', 4);
    if (iGroup == 1)
        hB = boxchart(iGroup*ones(sum(sampleGroup==iGroup), 1), yValues(sampleGroup==iGroup), 'BoxFaceColor', cmap(iGroup,:), 'LineWidth', 2, 'MarkerStyle', 'none', 'WhiskerLineColor', cmap(iGroup,:)); %  % , 'Color', cmap(iGroup,:)
    end
end

xAlt = xValues(iSampleMutated);
plot(xAlt, yValues(iSampleMutated), 'o', 'MarkerFaceColor', (1+sColours.mutated)/2, 'Color', sColours.mutated, 'MarkerSize', 6);

set(gca, 'XTick', []);
grid on;


yLimVal = get(gca, 'YLim'); yLimVal(1) = 0; yLimVal(2) = max([yLimVal(2); yValues(iSampleMutated)*1.1]);
ylim(yLimVal);

xAltRelative = (xAlt-xLimVal(1))/(xLimVal(2)-xLimVal(1));
yAltRelative = yValues(iSampleMutated)/yLimVal(2);

text(xLimVal(2), max(yValues(iSampleMutated)), 'MUT', 'Color', sColours.mutated, 'HorizontalAlignment', 'center', 'FontSize', fontSize-2);
text(xLimVal(2), median(yValues(sampleGroup==1), 'omitnan'), 'WT', 'Color', sColours.WT, 'HorizontalAlignment', 'center', 'FontSize', fontSize-2);

xlim(xLimVal);

ylabel(sprintf('{\\it%s} %s', geneName, yLabelText));
set(gca, 'FontSize', fontSize);
set(gca,'YColor',.5*[1,1,1]);

