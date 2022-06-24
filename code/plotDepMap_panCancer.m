function plotDepMap_panCancer(dataDepMap, sColours)

fontSize = 12;
fontSizeSmaller = fontSize - 2;
% figure
% hold on;
% yValues = dataDepMap.vDepMapGenes_pAbove50PerGene; yLabelText = 'Dependent cell lines (%)'; %saveName = 'pAbove50';
yValues = dataDepMap.vDepMapGenes_average; yLabelText = 'Average gene dependency'; %saveName = 'average';
groups = dataDepMap.groups;
labels = dataDepMap.lstGenesDepMap;

% xValues = groups+randn(length(groups), 1)/8;
xValues = myJitter(groups, yValues);



cmap = [sColours.WT; sColours.mutated];
for iGroup = 1:2
    if (iGroup == 2)
        plot(xValues(groups==iGroup), yValues(groups==iGroup), 'ok', 'MarkerFaceColor', (1+cmap(iGroup,:))/2, 'Color', cmap(iGroup,:), 'MarkerSize', 4);
    end
    boxchart(iGroup*ones(sum(groups==iGroup), 1), yValues(groups==iGroup), 'BoxFaceColor', cmap(iGroup,:), 'LineWidth', 2, 'MarkerStyle', '+', 'MarkerColor', .6*[1,1,1], 'WhiskerLineColor', cmap(iGroup,:)); %  % , 'Color', cmap(iGroup,:)
end

% boxchart(groups, yValues, 'BoxFaceColor', .5*[1,1,1], 'LineWidth', 2, 'MarkerStyle', '+', 'MarkerColor', .5*[1,1,1], 'WhiskerLineColor', .5*[1,1,1]); grid on;

isOK = groups==2;
plot(xValues(isOK), yValues(isOK), 'o', 'MarkerFaceColor', (1+sColours.mutated)/2, 'Color', sColours.mutated);

xShift = .05;
isOK_alignTop = ismember(labels, {'VPS28', 'IKBKB'});
text(xValues(isOK_alignTop) + xShift/2, yValues(isOK_alignTop), labels(isOK_alignTop), 'HorizontalAlignment', 'left', 'FontSize', fontSizeSmaller, 'FontAngle', 'italic', 'VerticalAlignment', 'top');

isOK1 = ~isOK_alignTop & isOK & yValues > 0.15 & xValues > 2;
text(xValues(isOK1) + xShift, yValues(isOK1), labels(isOK1), 'HorizontalAlignment', 'left', 'FontSize', fontSizeSmaller, 'FontAngle', 'italic');
isOK1 = isOK & yValues > 0.15 & xValues <= 2;
text(xValues(isOK1) - xShift, yValues(isOK1), labels(isOK1), 'HorizontalAlignment', 'right', 'FontSize', fontSizeSmaller, 'FontAngle', 'italic');
% set(gca, 'XTick', [], 'FontSize', 16);
ylabel(yLabelText);
p = ranksum(yValues(groups==1), yValues(groups==2), 'tail', 'both');
% title(sprintf('Solid cancers: %s', getPValueAsTextShort(p))); 


xlim([.6,2.5]); %ylim([0, 1.1]);

text(1,-0.05, {'Control', 'genes'}, 'HorizontalAlignment', 'center');
text(2,-0.05, {'Driver-upregulated', 'genes'}, 'HorizontalAlignment', 'center');

set(gca, 'XTick', [], 'FontSize', fontSize);
title(sprintf('DepMap\n\\fontsize{10}\\color[rgb]{0.5,0.5,0.5}{\\itp = %s}\n', getPValueAsTextTimes(p)));

drawnow;
