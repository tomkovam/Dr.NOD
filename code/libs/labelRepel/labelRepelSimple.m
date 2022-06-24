function [xValuesText, yValuesText] = labelRepelSimple(xValues, yValues, labels, fontSize, ignoreDist, maxForce, shiftDist, scaleBorderForce)

% ignoreDist = 0.05; maxForce = 0.01; shiftDist = 0.015;
if (~exist('ignoreDist', 'var'))
    ignoreDist = 0.05;
end

if (~exist('maxForce', 'var'))
    maxForce = 0.01;
end

if (~exist('shiftDist', 'var'))
    shiftDist = 0.015;
end

if (~exist('scaleBorderForce', 'var'))
    scaleBorderForce = 1;
end
seed = 2; %333
rng(seed);

nLabels = length(labels); %nMarkers = nLabels;
if (size(labels, 2) > 1) % convert to a column
    labels = labels';
end
if (size(xValues, 2) > 1) % convert to a column
    xValues = xValues';
end
if (size(yValues, 2) > 1) % convert to a column
    yValues = yValues';
end
xValuesOriginal = xValues;
yValuesOriginal = yValues;
xLimVal = get(gca, 'XLim'); xMax = max(xLimVal); xLimValNorm = xLimVal/xMax; %xLength = xLimVal(2)-xLimVal(1); 
yLimVal = get(gca, 'YLim'); yMax = max(yLimVal); yLimValNorm = yLimVal/yMax; %yLength = yLimVal(2)-yLimVal(1);  
posMark = [xValuesOriginal/xMax, yValuesOriginal/yMax];
posText  = [xValues/xMax, yValues/yMax];
posBoth = [posText; posMark]; nBoth = size(posBoth, 1);
xTickVal = get(gca, 'XTick')/xMax;
yTickVal = get(gca, 'YTick')/yMax;

%%
% gcaPos = get(gca, 'Position'); 
% distancesToMarkers = pdist2([xValues, yValues], [xValuesOriginal, yValuesOriginal], 'euclidean');
% distancesToTexts = pdist2([xValues, yValues], [xValues, yValues], 'euclidean');

% hTexts = text(posText(:,1)*xMax, posText(:, 2)*yMax, labels, 'HorizontalAlignment', 'center', 'FontSize', fontSize); drawnow; %pause(0.5);
for iStep = 1:50
    hasChanged = false;
    for iLabel = 1:nLabels
        force = [0,0];        
        lst = [1:iLabel-1, iLabel+1:nBoth];
        minDist = min(sqrt((posBoth(iLabel,1)-posBoth(lst,1)).^2 + (posBoth(iLabel,2)-posBoth(lst,2)).^2)); % pdist2(posBoth(iLabel,:), posBoth(lst,:), 'euclidean'));
        if (minDist < ignoreDist) % || posBoth(iLabel,1) < xLimValNorm(1) || posBoth(iLabel,2) < yLimValNorm(1))
            for jElement = lst
                vectorFromElement = posBoth(iLabel, :)-posBoth(jElement,:);
                distanceFromElement = sqrt(sum(vectorFromElement.^2));
                if (distanceFromElement < 0.5)
                    if (sum(distanceFromElement) == 0)
                        vectorFromElement = rand(1,2)*ignoreDist/3;
                        distanceFromElement = sqrt(sum(vectorFromElement.^2));
                    end
                    force = force + (maxForce*exp(-distanceFromElement*30))*vectorFromElement/distanceFromElement; %force + (maxForce*(1-distanceFromElement))*vectorFromElement/distanceFromElement;
                end
            end
            for iBorder = 1:4
                if (iBorder == 1)
                    vectorFromElement = [posBoth(iLabel, 1) - xLimValNorm(1), 0];
                elseif (iBorder == 2)
                    vectorFromElement = [posBoth(iLabel, 1) - xLimValNorm(2), 0];
                elseif (iBorder == 3)
                    vectorFromElement = [0, posBoth(iLabel, 2) - yLimValNorm(1)];
                elseif (iBorder == 4)
                    vectorFromElement = [0, posBoth(iLabel, 2) - yLimValNorm(2)];
                end
                distanceFromElement = sqrt(sum(vectorFromElement.^2));
                if (distanceFromElement < 2*ignoreDist)
                    if (sum(distanceFromElement) == 0)
                        vectorFromElement = rand(1,2)*ignoreDist/3;
                        distanceFromElement = sqrt(sum(vectorFromElement.^2));
                    end
                    force = force + (scaleBorderForce*maxForce*exp(-distanceFromElement*30))*vectorFromElement/distanceFromElement;
                end
            end
            %             posBoth(iLabel,1) = max([xLimValNorm(1), min([xLimValNorm(2), posBoth(iLabel,1) + force(1)])]);
            %             posBoth(iLabel,2) = max([yLimValNorm(1), min([yLimValNorm(2), posBoth(iLabel,2) + force(2)])]);
            posBoth(iLabel,1) = posBoth(iLabel,1) + force(1);
            posBoth(iLabel,2) = posBoth(iLabel,2) + force(2);
            if (scaleBorderForce > 0)
                posBoth(iLabel,1) = max([xLimValNorm(1), min([xLimValNorm(2), posBoth(iLabel,1)])]);
                posBoth(iLabel,2) = max([yLimValNorm(1), min([yLimValNorm(2), posBoth(iLabel,2)])]);
            end
        end
        hasChanged = true;
        if (scaleBorderForce == 0) %% posBoth(iLabel,:)
            if (posBoth(iLabel,1) < xLimValNorm(1) || posBoth(iLabel,2) < yLimValNorm(1)) % outside axes   % || posBoth(iLabel,1) > xLimValNorm(2)  || posBoth(iLabel,2) > yLimValNorm(2)
                force = [0,0];
                for iXY = 1:(length(xTickVal)+length(yTickVal))
                    if (iXY <= length(xTickVal) && (posBoth(iLabel,2) < yLimValNorm(1)))    % below axes
                        iTick = iXY;                    vectorFromElement = [posBoth(iLabel, 1) - xTickVal(iTick), 0]; %distanceToPlot = posBoth(iLabel,2) - yLimValNorm(1);
                    elseif (iXY > length(xTickVal) && posBoth(iLabel,1) < xLimValNorm(1))   % left from axes
                        iTick = iXY-length(xTickVal);   vectorFromElement = [0, posBoth(iLabel, 2) - yTickVal(iTick)]; %distanceToPlot = posBoth(iLabel,1) - xLimValNorm(1);
                    end
                    distanceFromElement = sqrt(sum(vectorFromElement.^2));
                    if (distanceFromElement < ignoreDist)
                        if (sum(distanceFromElement) == 0)
                            vectorFromElement = rand(1,2)*ignoreDist/3;
                            distanceFromElement = sqrt(sum(vectorFromElement.^2));
                        end
                        %vectorFromElement(vectorFromElement == 0) = 2*vectorToPlot; % To move it back into the axes.
                        force = force + (maxForce*exp(-distanceFromElement*30))*vectorFromElement/distanceFromElement;
                    end
                end
                if ((posBoth(iLabel,2) < yLimValNorm(1)))    % below axes
                    vectorToPlot = [0, yLimValNorm(1) - posBoth(iLabel,2)];
                elseif (posBoth(iLabel,1) < xLimValNorm(1))   % left from axes
                    vectorToPlot = [xLimValNorm(1) - posBoth(iLabel,1), 0];
                end
                distanceFromElement = sqrt(sum(vectorToPlot.^2));
                force = force + 0.2*(2^distanceFromElement - 1)*vectorToPlot/distanceFromElement;
                        
                posBoth(iLabel,1) = posBoth(iLabel,1) + force(1);
                posBoth(iLabel,2) = posBoth(iLabel,2) + force(2);
                hasChanged = true;
            end
        end
    end
    if (exist('hTexts', 'var'))
        delete(hTexts);
    end
    hTexts = text(posBoth(1:nLabels,1)*xMax, posBoth(1:nLabels, 2)*yMax, labels, 'HorizontalAlignment', 'center', 'FontSize', fontSize); drawnow; %pause(0.5);
    if (~hasChanged)
        fprintf('Nothing has changed.\n');
        break
    end
end
%%
lineStarts = posBoth(1:nLabels,:);
lineEnds = posMark;
lineVectors = lineEnds-lineStarts;
lineLengths = sqrt(sum(lineVectors.^2, 2));

shiftLengths = lineLengths/2; shiftLengths(shiftLengths>shiftDist) = shiftDist; % min from lineLengths/2 and shiftDist
lineStartsNew = lineStarts + lineVectors.*repmat(shiftLengths./lineLengths, 1, 2);
lineEndsNew = lineEnds - lineVectors.*repmat(0.5*shiftLengths./lineLengths, 1, 2);
% lineLengthsNew = sqrt(sum((lineEndsNew-lineStartsNew).^2, 2));
% lineStartsNew = lineStartsNew(lineLengthsNew>shiftDist,:);
% lineEndsNew = lineEndsNew(lineLengthsNew>shiftDist,:);

if (exist('hLines', 'var'))
    delete(hLines);
end
hLines = plot([lineStartsNew(:,1)*xMax, lineEndsNew(:,1)*xMax]', [lineStartsNew(:,2)*yMax, lineEndsNew(:,2)*yMax]', '-', 'Color', 0.65*[1,1,1]);
xlim(xLimVal); 
ylim(yLimVal);
if (scaleBorderForce == 0)
    set(gca,'clipping','off')% turn off clippings
end
if (exist('hTexts', 'var'))
    delete(hTexts);
end
xValuesText = posBoth(1:nLabels,1)*xMax;
yValuesText = posBoth(1:nLabels,2)*yMax;


% hLines = plot([posBoth(1:nLabels,1)*xMax, posMark(:,1)*xMax]', [posBoth(1:nLabels,2)*yMax, posMark(:,2)*yMax]', '-');
%for jElement = 1:nMarkers
% vectorFromElement = posText(iLabel,:)-posMark(jElement, :);
% if (vectorFromElement == [0,0])
%     vectorFromElement = rand(1,2)/100;
% end
% forceSize = sqrt(sum(vectorFromElement.^2));
% force = force + vectorFromElement./(1+forceSize^3);
% end