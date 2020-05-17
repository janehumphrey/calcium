%%
% NAME: PREPARE AXES
% AUTHOR: JANE HUMPHREY (janehumphrey@outlook.com)

function [ax] = prepareAxes(style,xLabel,yLabel,xLim,yLim,zLim)

if nargin<1
    error('Not enough input arguments.');
end

hold on;
ax = gca;
if strcmp(style,'image')==1
    colormap(gray);
    axis off;
    axis equal;
    ax.YDir = 'reverse';
    if nargin>3&&~isempty(xLim)
        ax.XLim = xLim;
    end
    if nargin>4
        ax.YLim = yLim;
    end
elseif strcmp(style,'plot')==1
    ax.Box = 'off';
    if nargin>1&&~isempty(xLabel)
        ax.XLabel.String = xLabel;
    end
    if nargin>2&&~isempty(yLabel)
        ax.YLabel.String = yLabel;
    end
    if nargin>3&&~isempty(xLim)
        ax.XLim = xLim;
    end
    if nargin>4
        ax.YLim = yLim;
    end
    ax.LineWidth = 1.25;
    ax.FontSize = 15;
    ax.XLabel.FontSize = 20;
    ax.YLabel.FontSize = 20;
elseif strcmp(style,'surface')==1
    colormap(gray);
    axis off;
    ax.Position = ax.OuterPosition;
    if nargin>3&&~isempty(xLim)
        ax.XLim = xLim;
    end
    if nargin>4&&~isempty(yLim)
        ax.YLim = yLim;
    end
    if nargin>5
        ax.ZLim = zLim;
        ax.CLim = zLim;
    end
else
    error('Unknown figure style.');
end