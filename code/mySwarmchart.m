function h = mySwarmchart(xValues, yValues, MarkerEdgeColor, MarkerFaceColor)

if (exist('swarmchart', 'file'))
    h = swarmchart(xValues, yValues, 'MarkerEdgeColor', MarkerEdgeColor, 'MarkerFaceColor', MarkerFaceColor);
else
    h = swarmchart(xValues+rand(length(xValues), 1)/4-1/8, yValues, 'MarkerEdgeColor', MarkerEdgeColor, 'MarkerFaceColor', MarkerFaceColor);
end
