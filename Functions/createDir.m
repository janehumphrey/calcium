%%
% NAME: CREATE DIRECTORY
% AUTHOR: JANE HUMPHREY (janehumphrey@outlook.com)

function [resultsDir] = createDir(resultsDir,inputDir,codeUsed)

if nargin<1
    error('Not enough input arguments');
end

if isempty(resultsDir)
    if nargin<3
        error('Not enough input arguments');
    end
    dirNames = strsplit(inputDir,filesep);
    dataStr = '';
    for i = 4:length(dirNames)
        dataStr = strcat(dataStr,filesep,dirNames{i});
    end
    resultsDir = ['C:\Dropbox (Personal)\RESULTS\',codeUsed,dataStr];
end

firstDir = resultsDir;
dirCount = 0;
while isfolder(resultsDir)
    dirCount = dirCount+1;
    resultsDir = strcat(firstDir,'_',num2str(dirCount));
end
mkdir(resultsDir);