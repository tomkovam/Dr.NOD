function yAltRelative = plotMotif(matFrequency, fontSize, lineWidth, sequence_ref, sequence_alt, motifName, iPosMutation)


sColours.A = [18,151,72]/256; %[0,1,0];
sColours.C = [32,98,156]/256; %[0,0,1];
sColours.G = [250,177,48]/256; %[16,151,72]/256; %[1,1,0];
sColours.T = [216,39,59]/256; %[1,0,0];
% lineWidth = 5; 
markerSize = 1;

lstBases = {'A', 'C', 'G', 'T'};
xEpsilon = 0.05;
% yEpsilon = 0.005;
yEpsilon = 0.02;

hold on;

maxY = 1.65;
x1 = iPosMutation + 0.5*[-1,1];
y_bottom = 0*[1,1];
y_top = maxY*[1,1];
h = patch([x1, fliplr(x1)], [y_bottom, fliplr(y_top)], 'k', 'EdgeColor', 'none');  set(h, 'FaceAlpha', 0.1);

xText = 0.2;

nPositions = size(matFrequency, 2);
for iPosition = 1:nPositions
    %matFrequency(:,iPosition)
    [heightRows, orderRows] = sort(matFrequency(:,iPosition));

    %     uncertaintyH = -sum(heightRows .* log(heightRows));
    %     %e_n = 1/log(2) * (4-1)/2*n; % We don't have the value of n...
    %     informationContentR = log2(4)- (uncertaintyH);

    y1 = 0;
    for iRow = 1:4
        heightBase = heightRows(iRow);
        %heightBase = heightBase .* informationContentR;
        if (heightBase > yEpsilon)
            plotOneBase(orderRows(iRow), iPosition - 0.5, iPosition + 0.5, y1, y1 + heightBase, xEpsilon, yEpsilon, sColours.(lstBases{orderRows(iRow)}), lineWidth, markerSize);  
        end
        y1 = y1 + heightBase;
    end
end

if (length(motifName)>9)
    fontSizeLabel = fontSize-2;
else
    fontSizeLabel = fontSize;
end

text(xText, y1/2, strrep(motifName, '_', '\_'), 'Rotation', 90, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', fontSizeLabel); % {strrep(motifName, '_', '\_'), 'motif'}

% xLimValues = [0.5, nPositions + 0.5];
% % plot(xLimValues, (1+yEpsilon)*[1,1], '-w', 'LineWidth', lineWidth);
% h = patch([xLimValues, fliplr(xLimValues)], [(y1+5*yEpsilon)*[1,1], fliplr((y1+10*yEpsilon)*[1,1])], 'w', 'EdgeColor', 'w');  %  set(h, 'FaceAlpha', 0.1);

y1 = y1 + 0.1;
for iPosition = 1:nPositions
    heightBase = 0.2;
    base = sequence_alt(iPosition);
    plotOneBase(find(strcmp(lstBases, base)), iPosition - 0.5, iPosition + 0.5, y1, y1 + heightBase, xEpsilon, yEpsilon, sColours.(base), lineWidth, markerSize);  
end
yAlt = y1+heightBase/2;
text(xText, yAlt, 'alt', 'HorizontalAlignment', 'right');
y1 = y1 + heightBase + 0.1;

for iPosition = 1:nPositions
    heightBase = 0.2;
    base = sequence_ref(iPosition);
    plotOneBase(find(strcmp(lstBases, base)), iPosition - 0.5, iPosition + 0.5, y1, y1 + heightBase, xEpsilon, yEpsilon, sColours.(base), lineWidth, markerSize);  
end
text(xText, y1+heightBase/2, 'ref', 'HorizontalAlignment', 'right');
y1 = y1 + heightBase;

xLimValues = [0.4, nPositions + 0.6];
xlim(xLimValues); ylim([0,maxY]); 
% ylabel('Probability');
set(gca, 'XTick', 1:nPositions, 'TickLength', [0 0], 'FontSize', fontSize, 'XTickLabelRotation', 0,'XColor', 'none','YColor','none'); % , 'Clipping', 'off'
ax = gca;               % get the current axis
ax.Clipping = 'off';    % turn clipping off
drawnow;

yAltRelative = yAlt/maxY;
%%
    function plotOneBase(iBase, x1, x2, y1, y2, xEpsilon, yEpsilon, colour, lineWidth, markerSize)
        x1 = x1 + xEpsilon;
        x2 = x2 - xEpsilon;
        y1 = y1 + yEpsilon;
        y2 = y2 - yEpsilon;
        %lineWidthThick = lineWidth + 2*lineWidth*(y2-y1);
        lineWidthThick = lineWidth;
        if (iBase == 1) % A
            plot([x1,(x1+x2)/2], [y1,y2], '-', 'LineWidth', lineWidth, 'MarkerSize', markerSize, 'Color', colour);
            plot([x2,(x1+x2)/2], [y1,y2], '-', 'LineWidth', lineWidth, 'MarkerSize', markerSize, 'Color', colour);
            plot((x1+x2)/2, y2, 's', 'LineWidth', 1, 'MarkerSize', 2, 'Color', 'none', 'MarkerFaceColor', colour);
            plot([x1+(x2-x1)/4,x1+(x2-x1)*3/4], (y1+y2)/2*[1,1], '-', 'LineWidth', lineWidthThick, 'Color', colour);
        elseif (iBase == 4) % T
            plot((x1+x2)/2*[1,1], [y1,y2], '-', 'LineWidth', lineWidth, 'MarkerSize', markerSize, 'Color', colour);
            plot([x1,x2], y2*[1,1], '-', 'LineWidth', lineWidthThick, 'MarkerSize', markerSize, 'Color', colour);
        elseif (iBase == 2) % C
            p = 0.9;
            x2 = x2 + 1-p;
            rX = (x2 - x1)/2;
            rY = (y2 - y1)/2;
            xMiddle = x1 + rX;
            yBottom = y1;
            n = 1000;
            k = max([1, round(n*(1-p))]);


            th = linspace(-0, pi*2, n);
            ang = th(k:end-k+1);
            xp=rX*cos(ang);
            yp=rY*sin(ang);
            yp = yp - min(yp);
            plot(xMiddle+xp,yBottom+yp, '-', 'LineWidth', lineWidth, 'MarkerSize', markerSize, 'Color', colour);
        elseif (iBase == 3) % G
            p = 0.9;
            x2 = x2 + 1-p;
            rX = (x2 - x1)/2;
            rY = (y2 - y1)/2;
            xMiddle = x1 + rX;
            yBottom = y1;
            n = 1000;
            k = max([1, round(n*(1-p))]);


            th = linspace(-0, pi*2, n);
            ang = th(k:end-k+1);
            xp=rX*cos(ang);
            yp=rY*sin(ang);
            yp = yp - min(yp);
            plot(xMiddle+xp,yBottom+yp, '-', 'LineWidth', lineWidth, 'MarkerSize', markerSize, 'Color', colour);

            x3 = xMiddle + xp(end);
            y3 = yBottom + yp(end);


            plot(x3*[1,1], [y3,(y1+y2)/2], '-', 'LineWidth', lineWidth, 'MarkerSize', markerSize, 'Color', colour);
            plot([(x1+x2)/2,x3], ((y1+y2)/2)*[1,1], '-', 'LineWidth', lineWidthThick, 'MarkerSize', markerSize, 'Color', colour);
            plot(x3*[1,1], [y3,(y1+y2)/2], 'o', 'LineWidth', lineWidthThick, 'MarkerSize', 2, 'Color', 'none', 'MarkerFaceColor', colour);
            %plot([(x1+x2)/2,x3+xEpsilon*2], (.9*(y1+y2)/2)*[1,1], '-', 'LineWidth', lineWidthThick, 'MarkerSize', markerSize, 'Color', colour);
        end
    end
end