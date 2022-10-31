% function plotSupFigure_housekeeping_genes(tableGencodeGenes, tableTissuesWithPancancer, sResults)

tableHousekeepingGenes = readtable(sProperties.GENES_HOUSEKEEPING); % 'data/genes/Housekeeping_GenesHuman.csv' % https://housekeeping.unicamp.br/Housekeeping_GenesHuman.csv https://academic.oup.com/nar/article/49/D1/D947/5871367
lstHousekeeping = tableHousekeepingGenes.Gene_name;
tableGencodeGenes.isHousekeeping = ismember(tableGencodeGenes.geneSymbol,lstHousekeeping);

tableTissuesWithPancancerUsed = tableTissuesWithPancancer(:,{'iTissue', 'tissue', 'tissuePrint', 'nSamplesWGS', 'nSamplesWGSandRNA', 'pFisherCDG', 'enrichmentCDG', 'nExpected_candidateDrivers_CDGs', 'nObserved_candidateDrivers_CDGs'});

%%
nTissues = length(sResults);

emptyVector = NaN*tableTissuesWithPancancerUsed.enrichmentCDG;
tableTissuesWithPancancerUsed.expectedCandidateDriver = emptyVector;
tableTissuesWithPancancerUsed.observedCandidateDriver = emptyVector;
tableTissuesWithPancancerUsed.pCandidateDriver = emptyVector;
tableTissuesWithPancancerUsed.enrichmentCandidateDriver = emptyVector;
%
tableTissuesWithPancancerUsed.expectedCandidateHousekeeping = emptyVector;
tableTissuesWithPancancerUsed.observedCandidateHousekeeping = emptyVector;
tableTissuesWithPancancerUsed.pCandidateHousekeeping = emptyVector;
tableTissuesWithPancancerUsed.enrichmentCandidateHousekeeping = emptyVector;
%
tableTissuesWithPancancerUsed.expectedHypomutatedHousekeeping = emptyVector;
tableTissuesWithPancancerUsed.observedHypomutatedHousekeeping = emptyVector;
tableTissuesWithPancancerUsed.pHypomutatedHousekeeping = emptyVector;
tableTissuesWithPancancerUsed.enrichmentHypomutatedHousekeeping = emptyVector;

tableGencodeGenes.isExpressed = false & tableGencodeGenes.isUsedGene;
tableGencodeGenes.isHypomutated = false & tableGencodeGenes.isUsedGene;
tableGencodeGenes.isHypomutatedSolid = false & tableGencodeGenes.isUsedGene;

for iTissue = 1:nTissues
    pM = sResults{iTissue}.pM;
    %     pE = sResults{iTissue}.pE;
    sizeEffectM = sResults{iTissue}.sizeEffectM;
    %     sizeEffectE = sResults{iTissue}.sizeEffectE;
    %     qCombined = sResults{iTissue}.qCombined;
    %     isUP = sResults{iTissue}.isUP;
    isCandidate = sResults{iTissue}.isCandidate;
    isDriver = sResults{iTissue}.isDriver;
    isHypomutated = sizeEffectM<1 & pM<0.05;
    geneName = sResults{iTissue}.geneName;
    %
    isHousekeeping = ismember(geneName, lstHousekeeping);
    %
    [p, enrichment, nObserved, nExpected] = myFisherTest(isCandidate, isDriver, 'both', false);
    tableTissuesWithPancancerUsed.pCandidateDriver(iTissue) = p;
    tableTissuesWithPancancerUsed.enrichmentCandidateDriver(iTissue) = enrichment;
    tableTissuesWithPancancerUsed.observedCandidateDriver(iTissue) = nObserved;
    tableTissuesWithPancancerUsed.expectedCandidateDriver(iTissue) = nExpected;
%     tableTissuesWithPancancerUsed.pTextCandidateDriver{iTissue} = getPValueAsText(p);

    [p, enrichment, nObserved, nExpected] = myFisherTest(isCandidate, isHousekeeping, 'both', false); 
    tableTissuesWithPancancerUsed.pCandidateHousekeeping(iTissue) = p;
    tableTissuesWithPancancerUsed.enrichmentCandidateHousekeeping(iTissue) = enrichment;
    tableTissuesWithPancancerUsed.observedCandidateHousekeeping(iTissue) = nObserved;
    tableTissuesWithPancancerUsed.expectedCandidateHousekeeping(iTissue) = nExpected;

    [p, enrichment, nObserved, nExpected] = myFisherTest(isHypomutated, isHousekeeping, 'both', false);
    tableTissuesWithPancancerUsed.pHypomutatedHousekeeping(iTissue) = p;
    tableTissuesWithPancancerUsed.enrichmentHypomutatedHousekeeping(iTissue) = enrichment;
    tableTissuesWithPancancerUsed.observedHypomutatedHousekeeping(iTissue) = nObserved;
    tableTissuesWithPancancerUsed.expectedHypomutatedHousekeeping(iTissue) = nExpected;

    tableGencodeGenes.isExpressed(ismember(tableGencodeGenes.geneSymbol, geneName)) = true;
    tableGencodeGenes.isHypomutated(ismember(tableGencodeGenes.geneSymbol, geneName(isHypomutated))) = true;
    if (iTissue > 1)
        tableGencodeGenes.isHypomutatedSolid(ismember(tableGencodeGenes.geneSymbol, geneName(isHypomutated))) = true;
    end
end
%%
isUsedGene = tableGencodeGenes.isUsedGene;
isUsedSolid = tableGencodeGenes.isUsedSolid;
%%
fprintf('Candidate & driver:\n');
[p, enrichment, nObserved, nExpected] = myFisherTest(tableGencodeGenes.isCandidate(isUsedGene), tableGencodeGenes.isDriver(isUsedGene), 'both', false); 
iTissue = nTissues + 2;
tableTissuesWithPancancerUsed.pCandidateDriver(iTissue) = p;
tableTissuesWithPancancerUsed.enrichmentCandidateDriver(iTissue) = enrichment;
tableTissuesWithPancancerUsed.observedCandidateDriver(iTissue) = nObserved;
tableTissuesWithPancancerUsed.expectedCandidateDriver(iTissue) = nExpected;

[p, enrichment, nObserved, nExpected] = myFisherTest(tableGencodeGenes.isCandidateSolid(isUsedSolid), tableGencodeGenes.isDriver(isUsedSolid), 'both', false); 
iTissue = nTissues + 1;
tableTissuesWithPancancerUsed.pCandidateDriver(iTissue) = p;
tableTissuesWithPancancerUsed.enrichmentCandidateDriver(iTissue) = enrichment;
tableTissuesWithPancancerUsed.observedCandidateDriver(iTissue) = nObserved;
tableTissuesWithPancancerUsed.expectedCandidateDriver(iTissue) = nExpected;
%%
fprintf('Candidate & housekeeping:\n');
[p, enrichment, nObserved, nExpected] = myFisherTest(tableGencodeGenes.isCandidate(isUsedGene), tableGencodeGenes.isHousekeeping(isUsedGene), 'both', false);  % This could be used if needed: the candidate regulatory drivers are depleted in housekeeping genes in lymphomas.
iTissue = nTissues + 2;
tableTissuesWithPancancerUsed.pCandidateHousekeeping(iTissue) = p;
tableTissuesWithPancancerUsed.enrichmentCandidateHousekeeping(iTissue) = enrichment;
tableTissuesWithPancancerUsed.observedCandidateHousekeeping(iTissue) = nObserved;
tableTissuesWithPancancerUsed.expectedCandidateHousekeeping(iTissue) = nExpected;

[p, enrichment, nObserved, nExpected] = myFisherTest(tableGencodeGenes.isCandidateSolid(isUsedSolid), tableGencodeGenes.isHousekeeping(isUsedSolid), 'both', false);  %  In solid cancers, there is a non-significant enrichment of 1.2x though, driven by 3x brain, 3x breast, 1x ovary, 1x liver
iTissue = nTissues + 1;
tableTissuesWithPancancerUsed.pCandidateHousekeeping(iTissue) = p;
tableTissuesWithPancancerUsed.enrichmentCandidateHousekeeping(iTissue) = enrichment;
tableTissuesWithPancancerUsed.observedCandidateHousekeeping(iTissue) = nObserved;
tableTissuesWithPancancerUsed.expectedCandidateHousekeeping(iTissue) = nExpected;
%%
fprintf('Hypomutated & housekeeping:\n');
[p, enrichment, nObserved, nExpected] = myFisherTest(tableGencodeGenes.isHypomutated(isUsedGene), tableGencodeGenes.isHousekeeping(isUsedGene), 'both', false); 
iTissue = nTissues + 2;
tableTissuesWithPancancerUsed.pHypomutatedHousekeeping(iTissue) = p;
tableTissuesWithPancancerUsed.enrichmentHypomutatedHousekeeping(iTissue) = enrichment;
tableTissuesWithPancancerUsed.observedHypomutatedHousekeeping(iTissue) = nObserved;
tableTissuesWithPancancerUsed.expectedHypomutatedHousekeeping(iTissue) = nExpected;

[p, enrichment, nObserved, nExpected] = myFisherTest(tableGencodeGenes.isHypomutatedSolid(isUsedSolid), tableGencodeGenes.isHousekeeping(isUsedSolid), 'both', false); 
iTissue = nTissues + 1;
tableTissuesWithPancancerUsed.pHypomutatedHousekeeping(iTissue) = p;
tableTissuesWithPancancerUsed.enrichmentHypomutatedHousekeeping(iTissue) = enrichment;
tableTissuesWithPancancerUsed.observedHypomutatedHousekeeping(iTissue) = nObserved;
tableTissuesWithPancancerUsed.expectedHypomutatedHousekeeping(iTissue) = nExpected;
%%
tableGencodeGenes(tableGencodeGenes.isHypomutated,:)
%%
tableTissuesWithPancancerUsed = tableTissuesWithPancancerUsed([end,end-1,1:end-2],:);
nRows = size(tableTissuesWithPancancerUsed, 1);
%%
% fig = createMaximisedFigure(3, [0 0 30 30]); 
% fontSize = 12;
% nR = 3; nC = 2; iS = 1; xS = 0.85; yS = 0.8; xB = 0.1; yB = 0.05; xM = -0.03; yM = -0.01;
% 
% myGeneralSubplot(nR,nC,iS,xS,yS,xB,yB,xM,yM); hold on; iS = iS + 1; 
% plotObservedExpected(tableTissuesWithPancancerUsed, sColours, nRows, fontSize, 'CandidateDriver');
% 
% myGeneralSubplot(nR,nC,iS,xS,yS,xB,yB,xM,yM); hold on; iS = iS + 1; 
% plotObservedExpected(tableTissuesWithPancancerUsed, sColours, nRows, fontSize, 'CandidateHousekeeping');
% 
% myGeneralSubplot(nR,nC,iS,xS,yS,xB,yB,xM,yM); hold on; iS = iS + 1; 
% plotObservedExpected(tableTissuesWithPancancerUsed, sColours, nRows, fontSize, 'HypomutatedHousekeeping');


% fontSizeLetters = 26;
% dim = [.007 .99 .01 .01]; str = 'a'; annotation('textbox',dim,'String',str, 'FontSize', fontSizeLetters, 'EdgeColor','none', 'FontWeight','bold');
% dim = [.505 .99 .01 .01]; str = 'b'; annotation('textbox',dim,'String',str, 'FontSize', fontSizeLetters, 'EdgeColor','none', 'FontWeight','bold');
% dim = [.007 .58 .01 .01]; str = 'c'; annotation('textbox',dim,'String',str, 'FontSize', fontSizeLetters, 'EdgeColor','none', 'FontWeight','bold');
% dim = [.505 .58 .01 .01]; str = 'd'; annotation('textbox',dim,'String',str, 'FontSize', fontSizeLetters, 'EdgeColor','none', 'FontWeight','bold');
% dim = [.007 .23 .01 .01]; str = 'e'; annotation('textbox',dim,'String',str, 'FontSize', fontSizeLetters, 'EdgeColor','none', 'FontWeight','bold');
%%
fig = createMaximisedFigure(3, [0 0 30 10]); hold on;
plotObservedExpected(tableTissuesWithPancancerUsed, sColours, nRows, fontSize, 'CandidateHousekeeping', 'Housekeeping genes within targets');

mySaveAs(fig, imagesPath, 'SupFig_housekeeping', false, true);
savefig([imagesPath, 'SupFig_housekeeping.fig']);


function plotObservedExpected(tableTissuesWithPancancer, sColours, nRows, fontSize, analysisType, yLabelText)
        matValues = [tableTissuesWithPancancer.(['expected', analysisType]), tableTissuesWithPancancer.(['observed', analysisType])];
        xValues = (1:size(tableTissuesWithPancancer, 1))';
        yValues = max(matValues, [], 2);
        hB = bar(matValues, 'EdgeColor', 'flat', 'FaceColor', 'flat');
        hB(1).CData = sColours.gray;
        hB(2).CData = sColours.darkRed;
        text(xValues, 3 + yValues, strcat(num2str((tableTissuesWithPancancer.(['enrichment', analysisType])), '%.1fx')), 'HorizontalAlignment', 'center', 'FontSize', fontSize-2, 'Color', .5*[1,1,1]); % {'n = '} insetad of 'n = ' will keep the space in there!
        text(xValues, 5 + yValues, arrayfun(@getPValueStarsAsText, tableTissuesWithPancancer.(['p', analysisType]), 'UniformOutput', false), 'HorizontalAlignment', 'center', 'FontSize', fontSize-2, 'Color', .5*[1,1,1]); % {'n = '} insetad of 'n = ' will keep the space in there!
        
        maxVal = max(5 + yValues); yGap = maxVal/20;
        yVal = 0.2*yGap + maxVal;        ylim([0, yVal]);
        %         yVal = 1.5*yGap + maxVal;
        %         text(xValues+.3, yVal+0*xValues, num2str(tableTissuesWithPancancer.nSamplesWGSandRNA, '%d'), 'HorizontalAlignment', 'right', 'FontSize', fontSize-2, 'Color', .5*[1,1,1]); % {'n = '} insetad of 'n = ' will keep the space in there!
        %         text(0, yVal, 'WGS+RNA', 'HorizontalAlignment', 'right', 'FontSize', fontSize-2, 'Color', .5*[1,1,1]); % {'n = '} insetad of 'n = ' will keep the space in there!
        %
        %         yVal = 3*yGap + maxVal;
        %         text(xValues+.3, yVal+0*xValues, num2str(tableTissuesWithPancancer.nSamplesWGS, '%d'), 'HorizontalAlignment', 'right', 'FontSize', fontSize-2, 'Color', .5*[1,1,1]); % {'n = '} insetad of 'n = ' will keep the space in there!
        %         text(0, yVal, 'WGS', 'HorizontalAlignment', 'right', 'FontSize', fontSize-2, 'Color', .5*[1,1,1]); % {'n = '} insetad of 'n = ' will keep the space in there!


        set(gca, 'XTick', xValues, 'XTickLabel', strrep(tableTissuesWithPancancer.tissuePrint, 'wo Blood', 'Solid'), 'XTickLabelRotation', 45, 'FontSize', fontSize, 'TickLength', [0 0]);
        ylabel(yLabelText);
        legend({'Expected', 'Observed'}, 'Location', 'NorthEast', 'FontSize', fontSize-2); legend boxoff %title(strrep(yTestName, '_', ' '));
        box off; xlim([0, nRows+1]);
    end
% end