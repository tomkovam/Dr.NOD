% myPrintMemory

myPrintMemory_sWhos = whos;
myPrintMemory_cellWhos = struct2cell(myPrintMemory_sWhos)';
myPrintMemory_fields = fieldnames(myPrintMemory_sWhos); 

myPrintMemory_tableWhos = table();
for myPrintMemory_iField = 1:7%length(fields)
    myPrintMemory_field = myPrintMemory_fields{myPrintMemory_iField};
    if (ismember(myPrintMemory_field, {'global', 'sparse', 'complex'}))
        myPrintMemory_field = ['is_', myPrintMemory_field];
    end
    if (strcmp(myPrintMemory_field, 'size'))
        myPrintMemory_tableWhos.(myPrintMemory_field) = cellfun(@(x) sprintf('%dx%d', x(1), x(2)), myPrintMemory_cellWhos(:, 2), 'UniformOutput', false);
    else
        myPrintMemory_tableWhos.(myPrintMemory_field) = myPrintMemory_cellWhos(:, myPrintMemory_iField);
    end
end

myPrintMemory_tableWhos.B = cell2mat(myPrintMemory_tableWhos.bytes);
myPrintMemory_tableWhos.GiB = floor(myPrintMemory_tableWhos.B/(1024^3));
myPrintMemory_tableWhos.MiB = floor(myPrintMemory_tableWhos.B/(1024^2));
myPrintMemory_tableWhos.KiB = floor(myPrintMemory_tableWhos.B/(1024^1));

myPrintMemory_tableWhos = sortrows(myPrintMemory_tableWhos,'B','ascend');

if (sum(myPrintMemory_tableWhos.MiB>=1)>0)
    myPrintMemory_tableWhos(myPrintMemory_tableWhos.MiB>=1, [1, 2, 4, 6, 9, 10, 11, 8])
end

fprintf('Total: %s B\n', num2sepNumStr(sum(myPrintMemory_tableWhos.B)));

clear myPrintMemory_sWhos myPrintMemory_cellWhos myPrintMemory_fields myPrintMemory_iField myPrintMemory_field myPrintMemory_tableWhos
