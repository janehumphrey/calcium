%%
% NAME: SELECT FILES
% AUTHOR: JANE HUMPHREY (janehumphrey@outlook.com)

function[fileList,fileExt,nFiles] = selectFiles(inputDir,fileType,nameIncludes,nameExcludes)

if nargin<2
    error('Not enough input arguments.');
end

fileList = cell(1000,1);
fileCount = 0;

dirContents = dir(inputDir);
nItems = size(dirContents,1);

for iItem = 3:nItems
    if dirContents(iItem).isdir==1
        subdir = [inputDir,filesep,dirContents(iItem).name];
        subdirContents = dir(subdir);
        nFiles = size(subdirContents,1);
        for iFile = 1:nFiles
            fileCount = fileCount+1;
            fileList{fileCount} = [dirContents(iItem).name,filesep,subdirContents(iFile).name];
        end
    else
        fileCount = fileCount+1;
        fileList{fileCount} = dirContents(iItem).name;
    end
end

index = true(fileCount,1);
for iFile = 1:fileCount
    name = fileList{iFile};
    if ~contains(name,fileType)
        index(iFile) = false;
    elseif nargin>2&&~contains(name,nameIncludes)
        index(iFile) = false;
    elseif nargin>3&&contains(name,nameExcludes)
        index(iFile) = false;
    end
end

fileList = fileList(index);
nFiles = sum(index);
fileExt = cell(nFiles,1);
for iFile = 1:nFiles
    dots = strfind(fileList{iFile},'.');
    nameEnd = dots(1)-1;
    fileExt{iFile} = fileList{iFile}(nameEnd+1:end);
    fileList{iFile} = fileList{iFile}(1:nameEnd);
end

if nFiles==0
    error('No files selected.');
end