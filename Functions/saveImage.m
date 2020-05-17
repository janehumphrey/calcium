%%
% NAME: SAVE IMAGE
% AUTHOR: JANE HUMPHREY (janehumphrey@outlook.com)

function [] = saveImage(fileBase)

if nargin<1
    error('Not enough input arguments.');
end

figFile = [fileBase,'.fig'];
savefig(figFile);
pngFile = [fileBase,'.png'];
export_fig(pngFile,'-png');