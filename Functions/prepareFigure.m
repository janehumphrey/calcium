%%
% NAME: PREPARE FIGURE
% AUTHOR: JANE HUMPHREY (janehumphrey@outlook.com)

function [fig] = prepareFigure(aspectRatio,scale)

figure;
fig = gcf;
fig.Color = 'w';

oldXPos = fig.Position(1);
oldYPos = fig.Position(2);
oldWidth = fig.Position(3);
oldHeight = fig.Position(4);

if nargin<1||isempty(aspectRatio)
    newHeight = oldHeight;
else
    newHeight = oldWidth/aspectRatio;
end
if nargin<2
    finalWidth = oldWidth;
    finalHeight = newHeight;
else
    finalWidth = oldWidth*scale;
    finalHeight = newHeight*scale;
end
newXPos = oldXPos+(oldWidth-finalWidth)/2;
newYPos = oldYPos+oldHeight-finalHeight;

fig.Position(1) = newXPos;
fig.Position(2) = newYPos;
fig.Position(3) = finalWidth;
fig.Position(4) = finalHeight;