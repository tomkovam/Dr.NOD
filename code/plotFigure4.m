function plotFigure4(imagesPath, sColours, tableGencodeGenes, tableTissues_data1, dataDepMap, sProperties)

fig = createMaximisedFigure(4, [0 0 30 30]); 
nR = 5; nC = 4; iS = 1+nC; xS = 0.7; yS = 0.5; xB = 0.08; yB = 0.08; xM = -0.03; yM = -0.05;

myGeneralSubplot(nR,nC,iS,(nC-1.1)+xS,1+yS,xB,yB,xM,yM); iS = iS + nC; hold on;
plotBarsWithLiterature(sColours, tableGencodeGenes);

%%
yB = 0.1; yS = 0.65; 
for iExample = 1
    if (iExample == 1)
        tissueName = 'liver'; geneName = 'EZH2';
    elseif (iExample == 2)
        tissueName = 'liver'; geneName = 'STOML2';
    elseif (iExample == 3)
        tissueName = 'lung'; geneName = 'ALOXE3';
    end
        iTissue = find(strcmp(tableTissues_data1.tissue, tissueName));

        %myGeneralSubplot(nR,nC,iS,1+xS,yS,xB,yB,xM,yM); iS = iS + 2; axPos = get(gca, 'Position') % 0.0800    0.4800    0.4037    0.1235
        axes('Position', [0.08, 0.5635, 0.38, 0.04]); 
        x_mutatedUE_relative1 = plotGene_genomicView_zoomedOutAndIn(tissueName, tableTissues_data1.biosampleABC{iTissue}, geneName, sColours, false, sProperties); %plotGene_genomicView_zoomedOut(tissueName, tableTissues_data1.biosampleABC{iTissue}, geneName, sColours, false, true, tableGencodeGenes);
        axPos1 = get(gca, 'Position');

        axes('Position', [0.08, 0.48, 0.38, 0.04]); 
        [x_mutatedUE_relative2, colour] = plotGene_genomicView_zoomedOutAndIn(tissueName, tableTissues_data1.biosampleABC{iTissue}, geneName, sColours, true, sProperties); %plotGene_genomicView(tissueName, tableTissues_data1.biosampleABC{iTissue}, geneName, sColours, true, false);
        axPos2 = get(gca, 'Position');

        for iEnd = 1:2
            xa = [axPos1(1) + axPos1(3)*x_mutatedUE_relative1(iEnd), axPos2(1) + axPos2(3)*x_mutatedUE_relative2(iEnd)];
            ya = [axPos1(2), axPos2(2) + axPos2(4)];
            annotation('line',xa,ya, 'Color', colour, 'LineWidth', 1.5);
        end

        iS = iS + 2;

        myGeneralSubplot(nR,nC,iS,xS,yS,xB,yB,xM,yM); iS = iS + 1;
        plotGene_boxplot(tissueName, tableTissues_data1.biosampleABC{iTissue}, geneName, sColours, true);

        myGeneralSubplot(nR,nC,iS,xS,yS,xB,yB,xM,yM); iS = iS + 1;
        plotGene_survivalKaplanMeier(geneName, tissueName, sColours, sProperties);    
end

yB = 0.08; yS = 0.5; 
for iExample = 1:6
    if (iExample == 1)
        tissueName = 'colorectal'; geneName = 'IER3';
    elseif (iExample == 2)
        tissueName = 'colorectal'; geneName = 'ZNF521';
    elseif (iExample == 3)
        tissueName = 'lung'; geneName = 'ALOXE3';
    elseif (iExample == 4)
        tissueName = 'pancreas'; geneName = 'PARP2';
    elseif (iExample == 5)
        tissueName = 'lung'; geneName = 'PLAU';
    elseif (iExample == 6)
        tissueName = 'lung'; geneName = 'CLTC';
    end
    iTissue = find(strcmp(tableTissues_data1.tissue, tissueName));
    myGeneralSubplot(nR,nC,iS,xS,yS,xB,yB,xM,yM); iS = iS + 1;
    plotGene_boxplot(tissueName, tableTissues_data1.biosampleABC{iTissue}, geneName, sColours, true);
    if (iExample == 3)
        iS = iS + 1;
    end
end

myGeneralSubplot(nR,nC,iS,xS,.95+yS,xB,yB,xM,yM); hold on;
plotDepMap_panCancer(dataDepMap, sColours);


fontSizeLetters = 26;
dim = [.007 .99 .01 .01]; str = 'a'; annotation('textbox',dim,'String',str, 'FontSize', fontSizeLetters, 'EdgeColor','none', 'FontWeight','bold');
dim = [.007 .65 .01 .01]; str = 'b'; annotation('textbox',dim,'String',str, 'FontSize', fontSizeLetters, 'EdgeColor','none', 'FontWeight','bold');
dim = [.500 .65 .01 .01]; str = 'c'; annotation('textbox',dim,'String',str, 'FontSize', fontSizeLetters, 'EdgeColor','none', 'FontWeight','bold');
dim = [.740 .65 .01 .01]; str = 'd'; annotation('textbox',dim,'String',str, 'FontSize', fontSizeLetters, 'EdgeColor','none', 'FontWeight','bold');
dim = [.007 .42 .01 .01]; str = 'e'; annotation('textbox',dim,'String',str, 'FontSize', fontSizeLetters, 'EdgeColor','none', 'FontWeight','bold');
dim = [.740 .42 .01 .01]; str = 'f'; annotation('textbox',dim,'String',str, 'FontSize', fontSizeLetters, 'EdgeColor','none', 'FontWeight','bold');


mySaveAs(fig, imagesPath, 'Fig4', true, true);
savefig([imagesPath, 'Fig4.fig']);
