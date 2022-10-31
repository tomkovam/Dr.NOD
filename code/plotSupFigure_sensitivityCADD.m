function plotSupFigure_sensitivityCADD(imagesPath, sResPanCancerCrossCADD, lstMinCADD_PHRED, tableGencodeGenes)

nCDGs = sum(tableGencodeGenes.isDriver);
%%
fig = createMaximisedFigure(1); 
nR = 2; nC = 2; iS = 1; xS = 0.8; yS = 0.7; xB = 0.1; yB = 0.1; xM = -0.04; yM = -0.04;

lstTypes = 2:3; legendValues = {'Pan-cancer solid', 'Pan-cancer'};
xValues = lstMinCADD_PHRED;

for iType = 1:4
    myGeneralSubplot(nR,nC,iS,xS,yS,xB,yB,xM,yM); iS = iS + 1; hold on;
    if (iType == 1)
        yValues = sResPanCancerCrossCADD.nCandidateDrivers;
        yText = 'Number of NRD targets'; %Targets of non-coding regulatory drivers';
    elseif (iType == 2)
        %         yValues = sResPanCancerCrossCADD.enrichmentCDG;
        %         yText = 'Log_2 CDG enrichment';
        yValues = sResPanCancerCrossCADD.nObserved_candidateDrivers_CDGs;
        yText = 'NRD targets that are CDGs';
    elseif (iType == 3)
        yValues = 100*(sResPanCancerCrossCADD.nObserved_candidateDrivers_CDGs./sResPanCancerCrossCADD.nCandidateDrivers);
        yText = 'CDGs in NRD targets (%)'; % targets of non-coding regulatory drivers
    elseif (iType == 4)
        yValues = 100*(sResPanCancerCrossCADD.nObserved_candidateDrivers_CDGs./nCDGs);
        yText = ' NRD targets in all CDGs (%)';
    end
    plot(xValues, yValues(lstTypes,:)', 'o-'); %lsline;
    ylabel(yText); xlabel('CADD PHRED cut-off'); 
    legend(legendValues, 'Location', 'Best');
    set(gca, 'FontSize', 16); grid on; grid minor;
end


fontSizeLetters = 26;
dim = [.030 .99 .01 .01]; str = 'a'; annotation('textbox',dim,'String',str, 'FontSize', fontSizeLetters, 'EdgeColor','none', 'FontWeight','bold');
dim = [.510 .99 .01 .01]; str = 'b'; annotation('textbox',dim,'String',str, 'FontSize', fontSizeLetters, 'EdgeColor','none', 'FontWeight','bold');
dim = [.030 .52 .01 .01]; str = 'c'; annotation('textbox',dim,'String',str, 'FontSize', fontSizeLetters, 'EdgeColor','none', 'FontWeight','bold');
dim = [.510 .52 .01 .01]; str = 'd'; annotation('textbox',dim,'String',str, 'FontSize', fontSizeLetters, 'EdgeColor','none', 'FontWeight','bold');


mySaveAs(fig, imagesPath, 'SupFig_sensitivityAndPower.png', true, true);
savefig([imagesPath, 'SupFig_sensitivityAndPower.fig']);

