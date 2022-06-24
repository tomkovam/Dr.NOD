function plotFigure1(imagesPath, sColours, sResults, lstGenesStrongSupport)
%%

nTissues = length(sResults);

fig = createMaximisedFigure(1, [0 0 30 15]);
% nR = round(sqrt(nTissues)); nC = ceil(nTissues/nR);

nR = 2; nC = 4; iS = 1; xS = 0.75; yS = 0.75; xB = 0.05; yB = 0.1; xM = -0.02; yM = -0.02;

fontSize = 10;

for iTissue = 1:nTissues
    %subplot(nR,nC,iTissue);
    myGeneralSubplot(nR,nC,iTissue,xS,yS,xB,yB,xM,yM); hold on;
    
    showMoreLabels = false;
    if (mod(iTissue, nC)==1)
        doPlotYLabel = true;
    else
        doPlotYLabel = false;
    end
    if (iTissue>(nR-1)*nC)
        doPlotXLabel = true;
    else
        doPlotXLabel = false;
    end
    if (iTissue == nTissues)
        doPlotLegend = true;
    else
        doPlotLegend = false;
    end
    plotTissueScatter(sColours, sResults, iTissue, fontSize, showMoreLabels, doPlotXLabel, doPlotYLabel, doPlotLegend, lstGenesStrongSupport);
end
% mySaveAs(fig, imagesPath, ['fig1Ext_',tissueName, '_', biosampleABC, '.png'], false, true);
mySaveAs(fig, imagesPath, 'Fig1.png', true, true);
