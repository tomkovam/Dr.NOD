function [x1, y1] = mySmoothing(xValues, yValues, nSteps1, nSteps2, degree)

isOk = ~isnan(xValues) & ~isnan(yValues);
xValues = xValues(isOk);
yValues = yValues(isOk);
minVal = min(xValues);
xValues = xValues-minVal;

ts = timeseries(yValues, xValues);
x0 = unique(sort([xValues, linspace(xValues(1), xValues(end), nSteps1)]));
tsOut = resample(ts, x0); y0 = squeeze(tsOut.Data)';
x1 = linspace(xValues(1), xValues(end), nSteps2);
[p,~,mu] = polyfit(x0, y0, degree); 
y1 = polyval(p, x1, [], mu);
x1 = x1 + minVal;