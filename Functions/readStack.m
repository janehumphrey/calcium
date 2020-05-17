%%
% NAME: READ STACK
% AUTHOR: JANE HUMPHREY (janehumphrey@outlook.com)

function [fovRaw,width,height,nFrames] = readStack(file,nFrames)

if nargin<1
    error('Not enough input arguments.');
end

fileInfo = imfinfo(file);
width = fileInfo.Width;
height = fileInfo.Height;
if nargin<2||isnan(nFrames)
    nFrames = numel(fileInfo);
end

fovRaw = zeros(height,width,nFrames,'uint16');
for iFrame = 1:nFrames
    fovRaw(:,:,iFrame) = imread(file,iFrame,'Info',fileInfo);
end