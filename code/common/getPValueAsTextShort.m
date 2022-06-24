function pValueAsText = getPValueAsTextShort(pValue)

% if (nargin < 2)
%     pValCoef = 1;
% end
% 
% pValCorrected = pValue * pValCoef;

% if (pValue < 0.00001)
%     textPart = '< 0.00001';
% else
%     textPart = sprintf('%.5f', pValue);
% end

if (pValue == 0)
    textPart = '0';
elseif (pValue < 0.001)
    textPart = sprintf('%.0e', pValue);
else
    textPart = sprintf('%.1g', pValue);
end

nStars = 0;
if (pValue < 0.001)
    nStars = 3;
elseif (pValue < 0.01)
    nStars = 2;
elseif (pValue < 0.05)
    nStars = 1;
end
if (nStars > 0)
    pValueAsText = sprintf('%s (%s)', textPart, repmat('*', 1, nStars));
else
    pValueAsText = textPart;
end

