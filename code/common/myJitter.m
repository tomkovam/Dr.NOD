function xValuesJittered = myJitter(xValues, yValues)

nBinsY = 30;
maxValueX = .5;
% xValues = [1*ones(50, 1);2*ones(50, 1)];
% yValues = randn(100, 1);
%%
% [yValues, permValues] = sort(yValues);
% permValues = randperm(length(yValues));
% yValues = yValues(permValues);
% xValues = xValues(permValues);
% permValuesInverse = permValues;
% permValuesInverse(permValues) = (1:length(permValues))';
%%
xValuesUnique = unique(xValues(~isnan(xValues)));
xValuesJittered = xValues;
for xValue = xValuesUnique'
    isOK = xValues == xValue & ~isnan(yValues);
    xValuesJittered(isOK) = myJitterOneColumn(xValues(isOK), yValues(isOK), nBinsY, maxValueX);
end
%%
% xValuesJittered = xValuesJittered(permValuesInverse);

end
%%
function xValuesJittered = myJitterOneColumn(xValues, yValues, nBinsY, maxValueX)
edges = linspace(min(yValues), max(yValues), nBinsY+1);
yBinValues = discretize(yValues, edges);
binCounts = histcounts(yValues, edges);
%%
xStep = maxValueX / max(binCounts);
xJitterTmp = (0:xStep:maxValueX)';
xJitterTmp2 = (0:2*xStep:maxValueX)';
xValuesJittered = xValues;
for iBin = unique(yBinValues)'
    isOKBin = yBinValues == iBin;
    if (sum(isOKBin)<5)
        xValuesJittered(isOKBin) = xValues(isOKBin) + xJitterTmp2(1:sum(isOKBin)) - xJitterTmp2(round(sum(isOKBin)/2));
    else
        xValuesJittered(isOKBin) = xValues(isOKBin) + xJitterTmp(1:sum(isOKBin)) - xJitterTmp(round(sum(isOKBin)/2));
    end
end
%%
if (length(yValues)>=20)
    [~, perm] = sort(yValues, 'descend'); nTop = 20;
    lstElements = perm(1:4:nTop); xValuesJittered(lstElements) = xValues(lstElements) + maxValueX/2;
    lstElements = perm(2:4:nTop); xValuesJittered(lstElements) = xValues(lstElements) - maxValueX/6;
    lstElements = perm(3:4:nTop); xValuesJittered(lstElements) = xValues(lstElements) + maxValueX/6;
    lstElements = perm(4:4:nTop); xValuesJittered(lstElements) = xValues(lstElements) - maxValueX/2;
end
%%
end


%%
% fig = createMaximisedFigure(1); hold on;
% plot(xValuesJittered, yValues, 'o');
% % cmap = lines(nBinsY);
% % for iBin = 1:nBinsY
% %     isOK = yBinValues == iBin;
% %     plot(xValuesJittered(isOK), yValues(isOK), 'o', 'Color', cmap(iBin,:), 'MarkerFaceColor', cmap(iBin,:));
% % end