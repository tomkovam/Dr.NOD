function [xAltRelative, yAltRelative] = plotGene_boxplot_forLogos(tissueName, biosampleABC, geneName, sColours, iSampleMutated)

% Created in saveForOneGeneVisualisation.m
load(['save/oneGene_', tissueName, '_', biosampleABC, '_', geneName], 'gene_pM', 'gene_pE', 'gene_qCombined', 'expressionPerSample', 'nSamples', 'sampleGroup', 'sampleGroupInclWoExpression');

%%
colors.n = .5*[1,1,1];
colors.group = 0*[1,1,1];
%%
% fig = createMaximisedFigure(1, [0 0 15 15]); axes('Position', [.15, .2, .8, .7]); 
fontSize = 10;
fontSizeSmaller = fontSize - 4;
lstGroups = {'WT', 'MUT'};

xLimVal = 1 + [-0.4, 0.5];

hold on;
xValues = 1+randn(nSamples, 1)/16;
% xValues = sampleGroup+randn(nSamples, 1)/16;
% yValues = log2(1+expressionPerSample); yLabelText = 'log_2 (1+FPKM_{uq})';
yValues = expressionPerSample; yLabelText = 'FPKM_{uq}';

cmap = [sColours.WT; sColours.mutated];
for iGroup = 1%:2
    plot(xValues(sampleGroup==iGroup), yValues(sampleGroup==iGroup), 'ok', 'MarkerFaceColor', (1+cmap(iGroup,:))/2, 'Color', cmap(iGroup,:), 'MarkerSize', 4);
    if (iGroup == 1)
        hB = boxchart(iGroup*ones(sum(sampleGroup==iGroup), 1), yValues(sampleGroup==iGroup), 'BoxFaceColor', cmap(iGroup,:), 'LineWidth', 2, 'MarkerStyle', 'none', 'WhiskerLineColor', cmap(iGroup,:)); %  % , 'Color', cmap(iGroup,:)
    end
end

% plot(xLimVal, yValues(iSampleMutated)*[1,1], '-', 'LineWidth', 2, 'Color', sColours.mutated);
% xAlt = 1+randn(1, 1)/16;
xAlt = xValues(iSampleMutated);
plot(xAlt, yValues(iSampleMutated), 'o', 'MarkerFaceColor', (1+sColours.mutated)/2, 'Color', sColours.mutated, 'MarkerSize', 6);
% text(xAlt, yValues(iSampleMutated), sprintf('%d', iSampleMutated));
% plot(xValues(iSampleMutated), yValues(iSampleMutated), 'sk', 'LineWidth', 2, 'MarkerSize', 6);

set(gca, 'XTick', []);
grid on;


yLimVal = get(gca, 'YLim'); yLimVal(1) = 0; yLimVal(2) = max([yLimVal(2), yValues(iSampleMutated)*1.1]);
ylim(yLimVal);

xAltRelative = (xAlt-xLimVal(1))/(xLimVal(2)-xLimVal(1));
yAltRelative = yValues(iSampleMutated)/yLimVal(2);

% text(xLimVal(2), yLimVal(2), sprintf('%s\n%s', insertColourIntoText('MUT', sColours.mutated), insertColourIntoText('WT', sColours.WT)), 'HorizontalAlignment', 'right', 'FontSize', fontSize-2); % [insertColourIntoText('MUT', sColours.mutated), '\n', insertColourIntoText('WT', sColours.WT)]
text(xLimVal(2), yValues(iSampleMutated), 'MUT', 'Color', sColours.mutated, 'HorizontalAlignment', 'center', 'FontSize', fontSize-2);
text(xLimVal(2), median(yValues(sampleGroup==1), 'omitnan'), 'WT', 'Color', sColours.WT, 'HorizontalAlignment', 'center', 'FontSize', fontSize-2);

    
% yVal = -1*yLimVal(2)/100; %min(yLimVal);
% text(1, yVal, sprintf('%s', lstGroups{1}), 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', 'FontSize', fontSize, 'Color', colors.group);
% text(2, yVal, sprintf('%s', lstGroups{2}), 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', 'FontSize', fontSize, 'Color', colors.group);
% 
% yVal = yVal - yLimVal(2)/5;
% text(.3, yVal, sprintf('WGS\nWGS+RNA'), 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', 'FontSize', fontSizeSmaller, 'Color', colors.n);
% text(1, yVal, sprintf('n = %d\nn = %d', sum(sampleGroupInclWoExpression==1), sum(sampleGroup==1)), 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', 'FontSize', fontSizeSmaller, 'Color', colors.n);
% text(2, yVal, sprintf('n = %d\nn = %d', sum(sampleGroupInclWoExpression==2), sum(sampleGroup==2)), 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', 'FontSize', fontSizeSmaller, 'Color', colors.n);

% text(1, 2*yVal, sprintf('\nn = %d', sum(sampleGroup==1)), 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', 'FontSize', fontSize, 'Color', colors.n);
% text(2, 2*yVal, sprintf('\nn = %d', sum(sampleGroup==2)), 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', 'FontSize', fontSize, 'Color', colors.n);
xlim(xLimVal);

ylabel([geneName, ' ', yLabelText]);
set(gca, 'FontSize', fontSize);
set(gca,'YColor',.5*[1,1,1]);
% set(gca,'YAxisLocation','right');

% if (printTissue)
%     titleText = sprintf('%s in %s', geneName, tissueName);
% else
%     titleText = geneName;
% end
% title(sprintf('%s\n\\fontsize{8}\\color[rgb]{0.5,0.5,0.5}{\\itp_E = %s, qCombined = %s}', titleText, getPValueAsTextTimes(gene_pE), getPValueAsTextTimes(gene_qCombined)));
% title(sprintf('{\\it%s} in %s', geneName, tissueName));

% if (abs(pE - mdl.Coefficients.pValue(2))>0.1), error('DIFFERENT pE values: %f %f', pE, mdl.Coefficients.pValue(2)); end

% mySaveAs(fig, imagesPath, [tissueName, '_', geneName, '_', biosampleABC,'.png'], false);
