function plotGene_survivalKaplanMeier(geneName, tissueName, sColours, sProperties)

fontSize = 12;
tableSummary = readtable(sProperties.PROTEIN_ATLAS_EXAMPLES, 'sheet', 'summary'); % 'data/proteinAtlasExamples.xlsx'
tableSummary.sheetName = strcat(tableSummary.gene, '_', tableSummary.tissue);
%%
% fig = createMaximisedFigure(2);
% for iRow = 1:9
iRow = find(strcmp(tableSummary.sheetName, [geneName, '_', tissueName]));
sheetName = tableSummary.sheetName{iRow};
% cutoff = tableSummary.cutoff(iRow);
p = tableSummary.p(iRow);
nSamplesLow = tableSummary.nSamplesLow(iRow);
nSamplesHigh = tableSummary.nSamplesHigh(iRow);
isNegativePrognostic = tableSummary.isNegativePrognostic(iRow);
%%
tableSurvival = readtable(sProperties.PROTEIN_ATLAS_EXAMPLES, 'sheet', sheetName); % 'data/proteinAtlasExamples.xlsx'
tableSurvival.Properties.VariableNames = {'Sample', 'Description', 'FPKM'};
isOK = contains(tableSurvival.Description, 'dead') | contains(tableSurvival.Description, 'alive'); if (sum(~isOK)>0), error('Some samples have missing survival status.'); end
tableSurvival.isDead = contains(tableSurvival.Description, 'dead');
tableSurvival.survivalDays = str2double(regexp(tableSurvival.Description,'(?<=, )\d+(?= days$)','once','match')); % For example: {'55 years, male, asian, stage:iiia, alive, 399 days'} --> 399
nSamples = size(tableSurvival, 1);
if (nSamples ~= (nSamplesLow + nSamplesHigh)), error('Incorrect number of samples'); end
%%
FPKM = tableSurvival.FPKM;
indexSample = (1:length(FPKM))';
survivalTimes = tableSurvival.survivalDays/365.242199; % years
isAlive = ~tableSurvival.isDead;
lineType = '-';
plotCI = false;
inPercents = false;
plotCensored = true;
%%
% subplot(3,3,iRow);
hold on;
isHigh = indexSample <= nSamplesHigh;
plotKMCurve(survivalTimes(isHigh), isAlive(isHigh), sColours.highExpression, lineType, plotCI, inPercents, plotCensored);
isLow = indexSample > nSamplesHigh;
plotKMCurve(survivalTimes(isLow), isAlive(isLow), sColours.lowExpression, lineType, plotCI, inPercents, plotCensored);
%%
xLimVal = [0, max(survivalTimes)*1.05];
xlim(xLimVal);
ylim([0,1]);

if (isNegativePrognostic)
    text(.05*xLimVal(2), 0.1, sprintf('high=%d', sum(isHigh)), 'Color', sColours.highExpression, 'HorizontalAlignment','left');
    text(.95*xLimVal(2), 0.9, sprintf('low=%d', sum(isLow)), 'Color', sColours.lowExpression, 'HorizontalAlignment','right');
end
xlabel('Years');
ylabel('Survival');
set(gca, 'FontSize', fontSize);
%%
% resultsA = logrank([survivalTimes(isHigh), isAlive(isHigh)], [survivalTimes(isLow), isAlive(isLow)], 0.05, 1);
% title(sprintf('%s: %.1g %.1g\nlow=%d, high=%d', tableSummary.gene{iRow}, resultsA.p, p, sum(isLow), sum(isHigh)));
title(sprintf('{\\it%s} in %s\n\\fontsize{10}\\color[rgb]{0.5,0.5,0.5}{\\itp = %s}', geneName, tissueName, getPValueAsTextTimes(p)));
drawnow;
% end