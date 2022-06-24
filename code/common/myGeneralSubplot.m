function output = myGeneralSubplot(nSubplotsRows, nSubplotsColumns, iSubplot, widthScaleParam, heightScaleParam, xShiftParam, yShiftParam, xRightMarginParam, yTopMarginParam, plotAxesParam)

% nSubplotsRows = 1;
% myGeneralSubplot(nR2,nC2,iS2,xS2,yS2,xB2,yB2,xM2,yM2)
% myGeneralSubplot(nR,nC,iS,xS,yS,xB,yB,xM,yM);
% 
% fig = figure(1); clf(fig);  set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 30 20],'units','normalized','outerposition',[0 0 1 1]);
% nR = 5; nC = 5; iS = 1; xS = 0.8; yS = 0.8; xB = 0.05; yB = 0.05; xM = 0.01; yM = 0.01;
% myGeneralSubplot(nR,nC,iS,xS,yS,xB,yB,xM,yM); iS = iS + 1; hold on;

widthScale = 0.9; %0.8
heightScale = 0.9; %0.85
if (nargin > 3)
    widthScale = widthScaleParam;
    heightScale = heightScaleParam;
end
xShift = 0.06;
yShift = 0.05;
if (nargin > 5)
    xShift = xShiftParam;
    yShift = yShiftParam;
end
xRightMargin = 0;
yTopMargin = 0;
if (nargin > 7)
    xRightMargin = xRightMarginParam;
    yTopMargin = yTopMarginParam;
end
plotAxes = true;
if (nargin > 9)
    plotAxes = plotAxesParam;
end


totalWidth = 1 - xRightMargin - xShift;
totalHeight = 1 - yTopMargin - yShift;

width = totalWidth/nSubplotsColumns;
height = totalHeight/nSubplotsRows;

% xCoord = (mod(iSubplot-1,nSubplotsColumns))*width;
% yCoord = floor((iSubplot-1)/nSubplotsColumns)*height;

xCoord = (mod(iSubplot-1,nSubplotsColumns))*width;
yCoord = totalHeight - (floor((iSubplot-1)/nSubplotsColumns)+1)*height;

positionVector = [xShift+xCoord, yShift+yCoord, width*widthScale, height*heightScale];

if (plotAxes)
    h = axes('Position', positionVector);
else
    h = 0;
end

output = positionVector; %[positionVector, h];