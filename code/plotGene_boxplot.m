function plotGene_boxplot(tissueName, biosampleABC, geneName, sColours, printTissue, exclusionType)
%%
if (~exist('exclusionType', 'var') || strcmp(exclusionType, 'excludePOLE_MSI'))
    suffix = '';
else
    suffix = ['_', exclusionType];
end
%%
% Created in saveForOneGeneVisualisation.m
load(['save/oneGene/oneGene_', tissueName, '_', biosampleABC, '_', geneName, suffix], 'gene_pM', 'gene_qCombined', 'gene_pE', ...
    'expressionPerSample', 'nSamples', 'sampleGroup', 'sampleGroupInclWoExpression', 'CNVperSample');

%%
colors.n = .5*[1,1,1];
colors.group = 0*[1,1,1];
%%
% fig = createMaximisedFigure(1, [0 0 15 15]); axes('Position', [.15, .2, .8, .7]); 
fontSize = 12;
fontSizeSmaller = fontSize - 4;
lstGroups = {'WT', 'MUT'};

hold on;
xValues = sampleGroup+randn(nSamples, 1)/16;
% yValues = log2(1+expressionPerSample); yLabelText = 'log_2 (1+FPKM_{uq})';
yValues = expressionPerSample; yLabelText = 'FPKM_{uq}';

cmap = [sColours.WT; sColours.mutated];
for iGroup = 1:2
    plot(xValues(sampleGroup==iGroup), yValues(sampleGroup==iGroup), 'ok', 'MarkerFaceColor', (1+cmap(iGroup,:))/2, 'Color', cmap(iGroup,:), 'MarkerSize', 4);
    hB = boxchart(iGroup*ones(sum(sampleGroup==iGroup), 1), yValues(sampleGroup==iGroup), 'BoxFaceColor', cmap(iGroup,:), 'LineWidth', 2, 'MarkerStyle', 'none', 'WhiskerLineColor', cmap(iGroup,:)); %  % , 'Color', cmap(iGroup,:)
end

% hB = boxplot(yValues, sampleGroup, 'symbol', ''); 
% set(hB, 'LineWidth', 2); 
% h = findobj(gca,'Tag','Box');
% for j=1:length(h)
%     patch(get(h(j),'XData'),get(h(j),'YData'),cmap(j,:),'FaceAlpha',.5);
% end

set(gca, 'XTick', []);
grid on;


yLimVal = get(gca, 'YLim'); yLimVal(1) = 0; ylim(yLimVal);
yVal = -1*yLimVal(2)/100; %min(yLimVal);
text(1, yVal, sprintf('%s', lstGroups{1}), 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', 'FontSize', fontSize, 'Color', colors.group);
text(2, yVal, sprintf('%s', lstGroups{2}), 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', 'FontSize', fontSize, 'Color', colors.group);

yVal = yVal - yLimVal(2)/5;
text(.3, yVal, sprintf('WGS\nWGS+RNA'), 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', 'FontSize', fontSizeSmaller, 'Color', colors.n);
text(1, yVal, sprintf('n = %d\nn = %d', sum(sampleGroupInclWoExpression==1), sum(sampleGroup==1)), 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', 'FontSize', fontSizeSmaller, 'Color', colors.n);
text(2, yVal, sprintf('n = %d\nn = %d', sum(sampleGroupInclWoExpression==2), sum(sampleGroup==2)), 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', 'FontSize', fontSizeSmaller, 'Color', colors.n);

% text(1, 2*yVal, sprintf('\nn = %d', sum(sampleGroup==1)), 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', 'FontSize', fontSize, 'Color', colors.n);
% text(2, 2*yVal, sprintf('\nn = %d', sum(sampleGroup==2)), 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', 'FontSize', fontSize, 'Color', colors.n);
xlim([0,2]+.5);

% ylabel(yLabelText);
ylabel(sprintf('{\\it%s} %s', geneName, yLabelText));
% set(gca, 'FontSize', fontSize);
set(gca, 'FontSize', 10);
set(gca,'YColor',.5*[1,1,1]);

% mdl =  fitglm([sampleGroup, CNVperSample], expressionPerSample, 'linear', 'Distribution', 'poisson', 'DispersionFlag', true) % Poisson with overdispersed count variable, as in https://www.nature.com/articles/s41467-019-13929-1 (quasi-Poisson family GLM)
% tmp1 = mdl.Coefficients.Estimate(1) + mdl.Coefficients.Estimate(2)*sampleGroup + mdl.Coefficients.Estimate(3)*CNVperSample;
% tmp2 = log(expressionPerSample);
% figure; plot(tmp1, tmp2, 'o');


% title(sprintf('%s in %s\n{\\itp_M = %s, p_E = %s\nqCombined = %s}', geneName, tissueName, getPValueAsTextShort(gene_pM), getPValueAsTextShort(gene_pE), getPValueAsTextShort(gene_qCombined)));
% title(sprintf('%s in %s\n\\color[rgb]{0.5,0.5,0.5}{\\itqCombined = %s}', geneName, tissueName, getPValueAsTextShort(gene_qCombined)));
if (printTissue)
    titleText = sprintf('{\\it%s} in %s', geneName, tissueName);
else
    titleText = sprintf('{\\it%s}', geneName);
end
title(sprintf('%s\n\\fontsize{10}\\color[rgb]{0.5,0.5,0.5}{\\itp_E = %s, q = %s}', titleText, getPValueAsTextTimes(gene_pE), getPValueAsTextTimes(gene_qCombined)));
% title(sprintf('{\\it%s} in %s', geneName, tissueName));

% if (abs(pE - mdl.Coefficients.pValue(2))>0.1), error('DIFFERENT pE values: %f %f', pE, mdl.Coefficients.pValue(2)); end

% mySaveAs(fig, imagesPath, [tissueName, '_', geneName, '_', biosampleABC,'.png'], false);
%% OLD

% lstUniqueEnhancers = find(matUniqueEnhancersGenes(:,iGene)); 
% 
% % tableUniqueEnhancers(lstUniqueEnhancers,1:10)
% 
% isOK = ismember(tableMutations_candidateOneTissue.iUniqueEnhancer, lstUniqueEnhancers);
% tableMutations_candidateOneTissue = tableMutations_candidateOneTissue(isOK, :);
% 
% if (excludeLowCADD)
%     tableMutations_candidateOneTissue = tableMutations_candidateOneTissue(tableMutations_candidateOneTissue.isHighCADD, :); % only high CADD
% end
% 
% 
% tableMutations_candidateOneTissue.expression = expressionPerSample(tableMutations_candidateOneTissue.iSample);
% 
% % tableMutations
% 
% if (true)
%     
%     tableMutations_candidateOneTissue = tableMutations_candidateOneTissue(~tableMutations_candidateOneTissue.isExcluded,:); % OLD: tableMutations(ismember(tableMutations.iSample, find(tableSamples.isExcluded)),:) = [];
% end

% sampleGroup = ones(nSamples, 1);
% sampleGroup(unique(tableMutations_candidateOneTissue.iSample)) = 2;
