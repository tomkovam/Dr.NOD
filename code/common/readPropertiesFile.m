function sProperties = readPropertiesFile(fileName)

fileID = fopen(fileName);
cellProperties = textscan(fileID,'%s %s','Delimiter','=','CommentStyle','#');
fclose(fileID);
parametersNames = cellProperties{1,1};
parametersValues = cellProperties{1,2};

sProperties = struct();
for iParameter = 1:length(parametersNames)
    status = 0;
    if (length(parametersValues{iParameter})<10)
        [parameterAsNumber, status] = str2num(parametersValues{iParameter});
    end
    if (status == 1)
        sProperties.(parametersNames{iParameter}) = parameterAsNumber;
    elseif (strcmp(parametersValues{iParameter}, 'true'))
        sProperties.(parametersNames{iParameter}) = true;
    elseif (strcmp(parametersValues{iParameter}, 'false'))
        sProperties.(parametersNames{iParameter}) = false;
    else
        sProperties.(parametersNames{iParameter}) = parametersValues{iParameter};
    end
end
