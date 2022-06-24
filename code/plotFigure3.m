function plotFigure3(imagesPath, sColours, tableMutations_candidate, tableTissues_data1, sResults, sProperties)

%%
tissueName = 'bloodLymphoid'; tissuePrint = 'blood';
iTissue = find(strcmp(tableTissues_data1.tissue, tissueName));
%%
tableTrinucleotides = readtable(sProperties.TABLE_TRINUCLEOTIDES); % 'data/tableTriNucl96.txt'
nPatterns = size(tableTrinucleotides, 1);
tableSignaturesCOSMIC = readtable(sProperties.COSMIC_SIGNATURES_SBS); % 'data/COSMIC_v3.2_SBS_GRCh37.txt'
tableSignaturesCOSMIC.refBase = cellfun(@(x) x(3), tableSignaturesCOSMIC.Type, 'UniformOutput', false);             % C from 'G[C>G]G'
tableSignaturesCOSMIC.altBase = cellfun(@(x) x(5), tableSignaturesCOSMIC.Type, 'UniformOutput', false);             % G from 'G[C>G]G'
tableSignaturesCOSMIC.trinucleotide = cellfun(@(x) x([1,3,7]), tableSignaturesCOSMIC.Type, 'UniformOutput', false); % GCG from 'G[C>G]G'
tableSignaturesCOSMIC = sortrows(tableSignaturesCOSMIC,{'refBase', 'altBase', 'trinucleotide'},'ascend');
tableSignaturesCOSMIC.patternName = strcat(tableSignaturesCOSMIC.trinucleotide, '>', tableSignaturesCOSMIC.altBase);
if (~isequal(tableSignaturesCOSMIC.patternName, tableTrinucleotides.patternName)), error('WRONG ORDER'); end
tableSignaturesCOSMIC = tableSignaturesCOSMIC(:,2:end-4);
%%
fig = createMaximisedFigure(3, [0 0 30 30]); 
nR = 6; nC = 4; iS = nC+1; xS = 0.7; yS = 0.5; xB = 0.08; yB = 0.08; xM = -0.03; yM = -0.03;

yB = 0.1;
myGeneralSubplot(nR,nC,iS,1+xS,.95+yS,xB,yB,xM,yM); iS = 3; hold on; fontSize = 12;
showMoreLabels = true; doPlotXLabel = true; doPlotYLabel = true; doPlotLegend = false;
plotTissueScatter(sColours, sResults, iTissue, fontSize, showMoreLabels, doPlotXLabel, doPlotYLabel, doPlotLegend);
set(gca, 'YTick', 0:4:16);

myGeneralSubplot(nR,nC,iS,1+xS,yS,xB,yB,xM,yM); iS = iS + nC; hold on; 
isOK = tableMutations_candidate.iTissue == iTissue & ~tableMutations_candidate.isExcluded & tableMutations_candidate.isHighCADD;
vector = histcounts(tableMutations_candidate.iPattern(isOK), 1:nPatterns+1);
plotMutationalSignature(vector, tissuePrint, tableTrinucleotides);

myGeneralSubplot(nR,nC,iS,1+xS,yS,xB,yB,xM,yM); iS = iS + 2; hold on; 
vector = tableSignaturesCOSMIC.SBS84;
textToShow = 'SBS84';
plotMutationalSignature(vector, textToShow, tableTrinucleotides);


yB = 0.08; 

lstGenes = {'MYC', 'SGK1', 'RFTN1', 'PIM1'}; 
for iExample = 1:length(lstGenes)
    geneName = lstGenes{iExample};

    myGeneralSubplot(nR,nC,iS,2+xS,yS,xB,yB,xM,yM); iS = iS + 3;
    plotGene_genomicView(tissueName, tableTissues_data1.biosampleABC{iTissue}, geneName, sColours, true, false, sProperties);

    myGeneralSubplot(nR,nC,iS,xS,yS,xB,yB,xM,yM); iS = iS + 1;
    plotGene_boxplot(tissueName, tableTissues_data1.biosampleABC{iTissue}, geneName, sColours, false);
end


fontSizeLetters = 26;
dim = [.007 .99 .01 .01]; str = 'a'; annotation('textbox',dim,'String',str, 'FontSize', fontSizeLetters, 'EdgeColor','none', 'FontWeight','bold');
dim = [.507 .99 .01 .01]; str = 'b'; annotation('textbox',dim,'String',str, 'FontSize', fontSizeLetters, 'EdgeColor','none', 'FontWeight','bold');
dim = [.507 .85 .01 .01]; str = 'c'; annotation('textbox',dim,'String',str, 'FontSize', fontSizeLetters, 'EdgeColor','none', 'FontWeight','bold');
dim = [.007 .68 .01 .01]; str = 'd'; annotation('textbox',dim,'String',str, 'FontSize', fontSizeLetters, 'EdgeColor','none', 'FontWeight','bold');
dim = [.727 .68 .01 .01]; str = 'e'; annotation('textbox',dim,'String',str, 'FontSize', fontSizeLetters, 'EdgeColor','none', 'FontWeight','bold');


mySaveAs(fig, imagesPath, 'Fig3.png', false, true);
return
%%
