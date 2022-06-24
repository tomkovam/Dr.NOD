function [h, x, f, numersAtRisk, flow, fup] = plotKMCurve(survivalTimes, isAlive, colourValue, lineType, plotCI, inPercents, plotCensored, timesAtRisk)

% islogical(censoredValues)
coef = 1;
if (inPercents)
    coef = 100;
end

maxTime = max(survivalTimes);
[f, x, flow, fup] = ecdf(survivalTimes, 'censoring', isAlive, 'function', 'survivor');
x = [0; x; maxTime]; f = [f(1); f; f(end)]; 
h = stairs(x, coef*f, lineType, 'Color', colourValue, 'LineWidth', 2); 
flow = [flow(1); flow; flow(end)]; fup = [fup(1); fup; fup(end)];
if (plotCI)
    stairs(x, coef*flow, ':', 'Color', colourValue); 
    stairs(x, coef*fup, ':', 'Color', colourValue);
end

censoredPatients = survivalTimes(isAlive); nCensored = length(censoredPatients);
isAlive = zeros(nCensored, 1);
for iCensored = 1:nCensored
    currentTime = censoredPatients(iCensored);
    isThisInterval = [x;1] >= currentTime & [1;x] <= currentTime;
    index = find(isThisInterval, 1, 'first');
    isAlive(iCensored) = f(index-1);
    %fprintf('currentTime %d: %f\n', currentTime, censoredValues(iCensored));
end
if (plotCensored)
    plot(censoredPatients, coef*isAlive, '+', 'Color', colourValue, 'LineWidth', 1.5);
end

if (~exist('timesAtRisk', 'var'))
    timesAtRisk = [];
end

numersAtRisk = NaN*timesAtRisk;
for iTime = 1:length(timesAtRisk)
    numersAtRisk(iTime) = sum(survivalTimes >= timesAtRisk(iTime));
end
