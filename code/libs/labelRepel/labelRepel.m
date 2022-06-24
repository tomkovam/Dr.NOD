function labelRepel(xValues, yValues, xShift, yShift, labels, fontSize, backgroundColor, extraPaddingRelative)

% nData = 30;
% xValues = 1:nData;
% yValues = rand(nData, 1);
% labels = cellstr(num2str(xValues', 'N=%-d'));
% fontSize = 12;
% fig = figure(1); clf;
% plot(xValues, yValues, 'o', 'MarkerFaceColor', 'c');

if (~exist('backgroundColor', 'var'))
    backgroundColor = [1,1,1];
end
if (~exist('extraPaddingRelative', 'var'))
    extraPaddingRelative = 0.005; % in 0-1 units of image, how large should be extra padding around labels. This will be converted to pixels, using the larger of xdim and ydim.
end

xScale = 5; %3
yScale = 5; %3
maxSteps = 2*200;
allowOutside = true;
% minDistance = 300;

figCurrent = gcf;
gcaCurrent = gca;
nLabels = length(labels);
if (size(labels, 2) > 1)
    labels = labels';
end
if (size(xValues, 2) > 1)
    xValues = xValues';
end
if (size(yValues, 2) > 1)
    yValues = yValues';
end
isCurrentHold = ishold;

hold on;
%%
% set(hLegend, 'Visible', 'off'); set(hLegend, 'Visible', 'on');
imgStruct = getframe(gcaCurrent); img3DMatrix = imgStruct.cdata/255;  %(2:end-1, 2:end-1, :)
imgIsBackgroundColor = (img3DMatrix(:,:,1) == backgroundColor(1) & img3DMatrix(:,:,2) == backgroundColor(2) & img3DMatrix(:,:,3) == backgroundColor(3));
% imgIsBackgroundColor = (flipud(imgIsBackgroundColor))';
xLimVal = get(gcaCurrent, 'XLim'); yLimVal = get(gcaCurrent, 'YLim');
hRectangle = rectangle('Position', [xLimVal(1), yLimVal(1), xLimVal(2)-xLimVal(1), yLimVal(2)-yLimVal(1)], 'FaceColor', 'w', 'EdgeColor', 'none');

xMiddle = mean(xLimVal);
yMiddle = mean(yLimVal);
gca_xMin = xLimVal(1);
gca_xWidth = xLimVal(2) - xLimVal(1);
gca_yMin = yLimVal(1);
gca_yWidth = yLimVal(2) - yLimVal(1);

image_xWidth = size(imgIsBackgroundColor, 2);
image_yWidth = size(imgIsBackgroundColor, 1);

labelWidth = NaN*ones(nLabels, 1);
labelHeight = NaN*ones(nLabels, 1);
tableLabels = table(labels, labelWidth, labelHeight); clear labelWidth labelHeight
for iLabel = 1:nLabels
    % iLabel = 1;
    hText = text(xMiddle, yMiddle, labels{iLabel}, 'FontSize', fontSize); drawnow;
    tmp_imgStruct = getframe(gcaCurrent); tmp_img3DMatrix = tmp_imgStruct.cdata/255; tmp_img3DMatrix = tmp_img3DMatrix(2:end-1, 2:end-1, :);
    tmp_imgIsBackgroundColor = (tmp_img3DMatrix(:,:,1) == backgroundColor(1) & tmp_img3DMatrix(:,:,2) == backgroundColor(2) & tmp_img3DMatrix(:,:,3) == backgroundColor(3));
    
    xProjection = max(~tmp_imgIsBackgroundColor, [], 1);
    yProjection = max(~tmp_imgIsBackgroundColor, [], 2);
    % size(xProjection)
    % size(yProjection)
    xMin = find(xProjection, 1, 'first');
    xMax = find(xProjection, 1, 'last');
    yMin = find(yProjection, 1, 'first');
    yMax = find(yProjection, 1, 'last');
    
    tableLabels.labelWidth(iLabel) = round(xMax - xMin);
    tableLabels.labelHeight(iLabel) = round(yMax - yMin);
    
    delete(hText);
    
    
    % fig2 = figure(223); clf;
    % imagesc(imgIsBackgroundColor-tmp_imgIsBackgroundColor);
    % imagesc(tmp_imgIsBackgroundColor);
    % figure(figCurrent);
end
tableLabels.xValue = round(1+image_xWidth*(xValues+xShift-gca_xMin)/gca_xWidth);
tableLabels.yValue = round(1+image_yWidth*(yValues+yShift-gca_yMin)/gca_yWidth);
tableLabels.xValueOriginal = tableLabels.xValue;
tableLabels.yValueOriginal = tableLabels.yValue;

delete(hRectangle);

hTexts = text(xValues, yValues, labels, 'FontSize', fontSize, 'HorizontalAlignment', 'center');
%%
% matrixWorking = imgIsBackgroundColor;
% fig2 = figure(222); clf; colormap('gray'); imagesc(flipud(matrixWorking'));

tableLabels.xValue = tableLabels.xValueOriginal;
tableLabels.yValue = tableLabels.yValueOriginal;


% We make sure that imgIsBackgroundColor has zeros along all the edges
imgIsBackgroundColor(1,:) = 0;
imgIsBackgroundColor(:,1) = 0;
imgIsBackgroundColor(end,:) = 0;
imgIsBackgroundColor(:,end) = 0;


tableLabels.isOk = false(nLabels, 1);
% isCurrentBackground = imgIsBackgroundColor;
% save
for iStep = 1:maxSteps/10
    changeMade = false; % if there has been a change in forces - we can stop otherwise
    %     disp(iStep);
    tableLabels.xForce = zeros(nLabels, 1);
    tableLabels.yForce = zeros(nLabels, 1);
    
    for iLabel = 1:nLabels
        %         disp(iLabel);
        isCurrentBackground = imgIsBackgroundColor;
        
        for jLabel = [1:iLabel-1,iLabel+1:nLabels]
            %             disp(jLabel);
            xIndices = round(tableLabels.xValue(jLabel)+(-tableLabels.labelWidth(jLabel)/2:tableLabels.labelWidth(jLabel)/2));%(0:tableLabels.labelWidth(jLabel)));
            yIndices = round(tableLabels.yValue(jLabel)+(-tableLabels.labelHeight(jLabel)/2:tableLabels.labelHeight(jLabel)/2));
            
            xIndices(xIndices<1 | xIndices > image_xWidth) = [];
            yIndices(yIndices<1 | yIndices > image_yWidth) = [];
            isCurrentBackground(image_yWidth-yIndices+1,  xIndices ) = 0;
        end
        
        % isCurrentBackground is eroded with a structural element of size
        % proportional to extraPaddingRelative. E.g., if image is 200x200
        % and extraPaddingRelative is 0.02, we'll apply a disk erosion with
        % radius of 4 pixels.
        sizePadding = round(extraPaddingRelative * max(size(isCurrentBackground,1),size(isCurrentBackground,2)));
        if (sizePadding > 0)
            erosionElement = strel('disk', sizePadding);
            isCurrentBackground = imerode(isCurrentBackground, erosionElement);
        end
        %         fig2 = figure(224); clf; colormap('gray'); imagesc(isCurrentBackground);
        %         figure(figCurrent);
        
        xIndices = round(tableLabels.xValue(iLabel)+(-tableLabels.labelWidth(iLabel)/2:tableLabels.labelWidth(iLabel)/2));%(0:tableLabels.labelWidth(iLabel)));
        yIndices = round(tableLabels.yValue(iLabel)+(-tableLabels.labelHeight(iLabel)/2:tableLabels.labelHeight(iLabel)/2));
        
        % Adding forces to keep the label inside the image
        if (~allowOutside)
            if (sum(xIndices<1) > 0)
                tableLabels.xForce(iLabel) = tableLabels.xForce(iLabel) + xScale*3;
                changeMade = true;
            elseif (sum(xIndices > image_xWidth) > 0)
                tableLabels.xForce(iLabel) = tableLabels.xForce(iLabel) - xScale*3;
                changeMade = true;
            elseif (sum(yIndices<1) > 0)
                tableLabels.yForce(iLabel) = tableLabels.yForce(iLabel) + yScale*3;
                changeMade = true;
            elseif (sum(yIndices > image_yWidth) > 0)
                tableLabels.yForce(iLabel) = tableLabels.yForce(iLabel) - yScale*3;
                changeMade = true;
            end
        end
        
        
        xIndices(xIndices<1 | xIndices > image_xWidth) = [];
        yIndices(yIndices<1 | yIndices > image_yWidth) = [];
        % matrixWorking(xIndices, yIndices) = 0;
        % fig2 = figure(223); clf; colormap('gray'); imagesc(flipud(matrixWorking'));
        
        currentBox = isCurrentBackground(image_yWidth-yIndices+1,  xIndices );
        xProjection = sum(~currentBox, 1); xMean = sum((1:length(xProjection)).*xProjection)/sum(xProjection); xMid = length(xProjection)/2 + 0.5;
        yProjection = sum(~currentBox, 2); yMean = sum((1:length(yProjection))'.*yProjection)/sum(yProjection); yMid = length(yProjection)/2 + 0.5; % check 0.5
        
        if (~tableLabels.isOk(iLabel) && iStep > round(maxSteps/10))
            break;
            xDirection = tableLabels.xValueOriginal(iLabel)-tableLabels.xValue(iLabel);
            yDirection = tableLabels.yValueOriginal(iLabel)-tableLabels.yValue(iLabel);
            forceLength = sqrt(xDirection^2 + yDirection^2);
            if (forceLength > 0)
                changeMade = true;
                xDirection = xDirection/forceLength;
                yDirection = yDirection/forceLength;
                tableLabels.xForce(iLabel) = tableLabels.xForce(iLabel) + xScale*xDirection;
                tableLabels.yForce(iLabel) = tableLabels.yForce(iLabel) + yScale*yDirection;
                tableLabels.isOk(iLabel) = false;
            end
        end
        
        if (sum(~currentBox(:)) > 0)
            tableLabels.isOk(iLabel) = false;
            changeMade = true;
            
            xDirection = xMean-xMid;
            yDirection = yMean-yMid;
            forceLength = sqrt(xDirection^2 + yDirection^2);
            if (forceLength == 0) % if forces sum to nothing, we add a random perturbation
                xDirection = randi(5)-2.5;
                yDirection = randi(5)-2.5;
                forceLength = sqrt(xDirection^2 + yDirection^2);
            end
            xDirection = xDirection/forceLength;
            yDirection = yDirection/forceLength;
            tableLabels.xForce(iLabel) = tableLabels.xForce(iLabel) - xScale*xDirection;
            tableLabels.yForce(iLabel) = tableLabels.yForce(iLabel) - yScale*yDirection;
            if (iStep == 1)
                scaleBigBox = 3;
                scaleBigBoxLabels = 0.5;
                xIndices = round(tableLabels.xValue(iLabel)+(-scaleBigBox*tableLabels.labelWidth(iLabel)/2:scaleBigBox*tableLabels.labelWidth(iLabel)/2));
                yIndices = round(tableLabels.yValue(iLabel)+(-scaleBigBox*tableLabels.labelHeight(iLabel)/2:scaleBigBox*tableLabels.labelHeight(iLabel)/2));
                
                xIndices(xIndices<1 | xIndices > image_xWidth) = [];
                yIndices(yIndices<1 | yIndices > image_yWidth) = [];
                currentBigBox = isCurrentBackground(image_yWidth-yIndices+1,  xIndices );
                %currentBigBoxSmoothed = imgaussfilt(1+currentBigBox, 2, 'FilterSize', 2*round([tableLabels.labelWidth(iLabel)/2, tableLabels.labelHeight(iLabel)/2]/2)+1);
                %currentBigBoxSmoothed = medfilt2(1+currentBigBox, 2*round([tableLabels.labelHeight(iLabel)/2, tableLabels.labelWidth(iLabel)/2]/2)+1);
                currentBigBoxSmoothed = imfilter(1.2*currentBigBox, ones(scaleBigBoxLabels*2*round([tableLabels.labelHeight(iLabel), tableLabels.labelWidth(iLabel)]/2)+1) / 25);
                [rowsMaxima,columnsMaxima] = find(currentBigBoxSmoothed==max(currentBigBoxSmoothed(:)));
                distances = (tableLabels.xValue(iLabel)-(min(xIndices)+columnsMaxima)).^2 + (tableLabels.yValue(iLabel)-(min(yIndices)+rowsMaxima)).^2;
                [~, indexMinDistance] = min(distances);
                tableLabels.xForce(iLabel) = min(xIndices)+columnsMaxima(indexMinDistance) - tableLabels.xValue(iLabel);
                tableLabels.yForce(iLabel) = min(yIndices)+rowsMaxima(indexMinDistance) - tableLabels.yValue(iLabel);                
            elseif (iStep == round(maxSteps/10))
                break;
                tableLabels.xForce(iLabel) = image_xWidth*rand(1) - tableLabels.xValue(iLabel);
                tableLabels.yForce(iLabel) = image_yWidth*rand(1) - tableLabels.yValue(iLabel);
            end
            %             if (iStep > maxSteps/3)
            %                 xDirection = randi(1)-1/2;
            %                 yDirection = randi(1)-1/2;
            %                 forceLength = sqrt(xDirection^2 + yDirection^2);
            %                 xDirection = xDirection/forceLength;
            %                 yDirection = yDirection/forceLength;
            %                 tableLabels.xForce(iLabel) = tableLabels.xForce(iLabel) - 5*xScale*xDirection;
            %                 tableLabels.yForce(iLabel) = tableLabels.yForce(iLabel) - 5*yScale*yDirection;
            %             end
        elseif (~changeMade)
            tableLabels.isOk(iLabel) = true;
        end
    end
    
    tableLabels.xValue = tableLabels.xValue + tableLabels.xForce;
    tableLabels.yValue = tableLabels.yValue + tableLabels.yForce;
    
    tableLabels.gca_xValue = gca_xMin + gca_xWidth*(tableLabels.xValue-1)/image_xWidth;
    tableLabels.gca_yValue = gca_yMin + gca_yWidth*(tableLabels.yValue-1)/image_yWidth;
    
    figure(figCurrent);
    delete(hTexts);
    hTexts = text(tableLabels.gca_xValue, tableLabels.gca_yValue, labels, 'FontSize', fontSize, 'HorizontalAlignment', 'center');
    if (exist('hStep', 'var'))
        delete(hStep);
    end
    hStep = text(gca_xMin, gca_yMin, sprintf('Step = %d', iStep), 'VerticalAlignment', 'bottom');
    drawnow;
    %pause(0.1);
    %     if (iStep == 1)
    %         a = 1;
    %     end
    if (~changeMade)
        break;
    end
end

% for iLabel = find(tableLabels.isOk')
%
% end
%%
% tableLabels.distanceFromPoint = sqrt((tableLabels.xValue - tableLabels.gca_xValue).^2 + (tableLabels.yValue - tableLabels.gca_yValue).^2);
tableLabels.xDistance = abs(xValues - tableLabels.gca_xValue);
tableLabels.yDistance = abs(yValues - tableLabels.gca_yValue);
tableLabels.isVeryDistant = (tableLabels.xDistance > gca_xWidth/20 | tableLabels.yDistance > gca_yWidth/20);
tableLabels.isDistant = (tableLabels.xDistance > gca_xWidth/35 | tableLabels.yDistance > gca_yWidth/35);
delete(hTexts);
tableLabels.isOk(:) = true;
% tableLabels.isOk(~tableLabels.isVeryDistant) = true;

% linePercentageGCA = 50;

if (false)
    for iLabel = find(tableLabels.isDistant & tableLabels.isOk)'
        %     if (tableLabels.isDistant(iLabel)) %tableLabels.distanceFromPoint(iLabel) > minDistance)
        xDirection = tableLabels.gca_xValue(iLabel) - xValues(iLabel);
        yDirection = tableLabels.gca_yValue(iLabel) - yValues(iLabel);
        forceLength = sqrt(xDirection^2 + yDirection^2);
        %         xStart = xValues(iLabel) + (gca_xWidth/linePercentageGCA)*xDirection/forceLength;
        %         yStart = yValues(iLabel) + (gca_yWidth/linePercentageGCA)*yDirection/forceLength;
        %         xEnd = xValues(iLabel) + xDirection - (gca_xWidth/linePercentageGCA)*xDirection/forceLength;
        %         yEnd = yValues(iLabel) + yDirection - (gca_yWidth/linePercentageGCA)*yDirection/forceLength;
        xStart = xValues(iLabel) + 0.3*xDirection/forceLength;
        yStart = yValues(iLabel) + 0.3*yDirection/forceLength;
        xEnd = xValues(iLabel) + xDirection - 0.3*xDirection/forceLength;
        yEnd = yValues(iLabel) + yDirection - 0.3*yDirection/forceLength;
        plot([xStart, xEnd], [yStart, yEnd], '-', 'Color', 0.5*[1,1,1]);
        %     end
    end
end
hTexts = text(tableLabels.gca_xValue(tableLabels.isOk), tableLabels.gca_yValue(tableLabels.isOk), labels(tableLabels.isOk), 'FontSize', fontSize, 'HorizontalAlignment', 'center');

%%
figure(figCurrent);
% delete(hTexts);
delete(hStep);
if (~isCurrentHold)
    hold off;
end

% fig2 = figure(222);
% imagesc(imgIsBackgroundColor);
% figure(figCurrent);

