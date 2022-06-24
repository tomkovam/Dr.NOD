function createDir(myPath)

if (~exist(myPath,'dir'))
    disp(['Creating directory ',myPath,'...']);
    mkdir(myPath);
end