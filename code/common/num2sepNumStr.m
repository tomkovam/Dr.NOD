function result = num2sepNumStr(inputNumber)

if (isscalar(inputNumber))
    result = sprintf(',%c%c%c',fliplr(num2str(inputNumber)));
    result = fliplr(result(2:end));
else
    result = cell(length(inputNumber), 1);
    for iElement = 1:length(inputNumber)
        tmp = sprintf(',%c%c%c',fliplr(num2str(inputNumber(iElement))));
        tmp = fliplr(tmp(2:end));
        result{iElement} = tmp;
    end
end