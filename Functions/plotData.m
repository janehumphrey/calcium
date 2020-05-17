%%
% NAME: PLOT DATA
% AUTHOR: JANE HUMPHREY (janehumphrey@outlook.com)

function [data] = plotData(style,xValues,yValues,colour,lineStyle,lineWidth,markerStyle,markerFill,markerSize,...
    textArray,textSize,labels)

if nargin<3
    error('Not enough input arguments.');
end
if isempty(xValues)
    if strcmp(style,'strip')==1||strcmp(style,'jitter')==1
        xValues = repmat(1:size(yValues,2),size(yValues,1),1);
    else
        xValues = (1:size(yValues,1))';
    end 
end
if nargin<4||isempty(colour)
    if strcmp(style,'diff')==1||strcmp(style,'horz')==1||strcmp(style,'vert')==1
        colour = [0,0,0];
    else
        colour = [1,0.5,0];
    end
end
if nargin<5||isempty(lineStyle)
    lineStyle = '-';
end
if nargin<6||isempty(lineWidth)
    if strcmp(style,'line')==1
        lineWidth = 3;
    else
        lineWidth = 2;
    end
end
if nargin<7||isempty(markerStyle)
    markerStyle = 'o';
end
if nargin<8||isempty(markerFill)
    markerFill = 'white';
end
if nargin<9||isempty(markerSize)
    if strcmp(style,'diff')==1
        markerSize = 60;
    else
        markerSize = 30;
    end
end
if nargin<12||isempty(labels)
    labels = [];
end
if strcmp(style,'points')==1||strcmp(style,'strip')==1||strcmp(style,'jitter')==1
    if markerSize<50
        lineWidth = 1.5;
    end
end

ax = gca;
if strcmp(style,'line')==1
    data = plot(xValues,yValues,lineStyle,'Color',colour,'LineWidth',lineWidth);
elseif strcmp(style,'points')==1
    if strcmp(markerFill,'white')==1
        for iSet = 1:size(yValues,2)
            data = scatter(xValues(:,iSet),yValues(:,iSet),markerSize,colour,markerStyle,'LineWidth',lineWidth,...
                'MarkerFaceColor',[1,1,1]);
        end
    elseif strcmp(markerFill,'colour')==1
        for iSet = 1:size(yValues,2)
            data = scatter(xValues(:,iSet),yValues(:,iSet),markerSize,colour,'filled',markerStyle,'LineWidth',...
                lineWidth);
        end
    elseif strcmp(markerFill,'none')==1
        for iSet = 1:size(yValues,2)
            data = scatter(xValues(:,iSet),yValues(:,iSet),markerSize,colour,markerStyle,'LineWidth',lineWidth);
        end
    else
        error('Unknown marker fill.');
    end
elseif strcmp(style,'strip')==1
    for iSet = 1:size(yValues,2)
        data = scatter(xValues(:,iSet),yValues(:,iSet),markerSize,colour,markerStyle,'LineWidth',lineWidth,...
            'MarkerFaceColor',[1,1,1]);
    end
    ax.XTick = 1:size(yValues,2);
    ax.XTickLabel = labels;
    ax.XTickLabelRotation = 60;
elseif strcmp(style,'jitter')==1
    if size(yValues,1)>3
        jitter = 'on';
    else
        jitter = 'off';
    end
    for iSet = 1:size(yValues,2)
        data = scatter(xValues(:,iSet),yValues(:,iSet),markerSize,colour,markerStyle,'LineWidth',lineWidth,...
            'MarkerFaceColor',[1,1,1],'jitter',jitter);
    end
    ax.XTick = 1:size(yValues,2);
    ax.XTickLabel = labels;
    ax.XTickLabelRotation = 60;
elseif strcmp(style,'text')==1
    data = text(xValues,yValues,textArray,'Color',colour,'FontSize',textSize,'FontWeight','bold',...
        'HorizontalAlignment','center');
elseif strcmp(style,'hist')==1
    data = histogram(yValues,'BinMethod','sqrt','LineWidth',1.5,'EdgeColor',[0.15,0.15,0.15],'FaceColor',colour,...
        'FaceAlpha',1);
elseif strcmp(style,'box')==1
    if strcmp(ax.XLimMode,'manual')==1
        xManual = true;
        xLim = ax.XLim;
    else
        xManual = false;
    end
    if strcmp(ax.YLimMode,'manual')==1
        yManual = true;
        yLim = ax.YLim;
    else
        yManual = false;
    end
    boxplot(yValues,'Colors',[0,0,0]);
    boxHandle = ax.Children.Children;
    nLines = size(boxHandle,1);
    for iLine = 1:nLines
        boxHandle(iLine).LineWidth = lineWidth;
        if strcmp(boxHandle(iLine).Tag,'Box')==1
            data = fill(boxHandle(iLine).XData,boxHandle(iLine).YData,colour);
            boxHandle(iLine).LineStyle = 'none';
            data.LineStyle = 'none';
            uistack(data,'bottom');
        elseif strfind(boxHandle(iLine).Tag,'Whisker')>0
            boxHandle(iLine).LineStyle = ':';
        elseif strcmp(boxHandle(iLine).Tag,'Outliers')==1
            boxHandle(iLine).Marker = 'x';
            boxHandle(iLine).MarkerSize = 10;
            boxHandle(iLine).MarkerEdgeColor = [0,0,0];
        end
    end
    ax.Box = 'off';
    if xManual==true
        ax.XLim = xLim;
    end
    if yManual==true
        ax.YLim = yLim;
    end
    ax.XTick = 1:size(yValues,2);
    ax.XTickLabel = labels;
    ax.XTickLabelRotation = 60;
elseif strcmp(style,'bar')==1
    data = bar(xValues,yValues,'FaceColor','flat','CData',colour,'LineWidth',1.5,'EdgeColor',[0.15,0.15,0.15]);
    ax.XTick = xValues;
    ax.XTickLabel = labels;
    ax.XTickLabelRotation = 60;
    data.BaseLine.LineWidth = 1.5;
    uistack(data,'bottom');
elseif strcmp(style,'error')==1
    errorbar(xValues,yValues(:,1),yValues(:,2),'LineStyle','none','Color',[0,0,0],'LineWidth',lineWidth,'CapSize',8);
elseif strcmp(style,'diff')==1
    if ~isvector(yValues)
        error('Data isn''t a vector.');
    end
    for iSet = 1:size(yValues,1)
        data = scatter(xValues(iSet),yValues(iSet),markerSize,colour,'filled',markerStyle,'LineWidth',lineWidth);
    end
elseif strcmp(style,'horz')==1
    if ~isvector(yValues)
        error('Data isn''t a vector.');
    end
    for iSet = 1:size(yValues,1)
        data = plot([xValues(iSet)-0.25,xValues(iSet)+0.25],[yValues(iSet),yValues(iSet)],lineStyle,'Color',colour,...
            'LineWidth',lineWidth);
    end
elseif strcmp(style,'vert')==1
    if size(yValues,2)~=2
        error('Data isn''t a 2-column matrix.');
    end
    for iSet = 1:size(yValues,1)
        data = plot([xValues(iSet),xValues(iSet)],[yValues(iSet,1),yValues(iSet,2)],lineStyle,'Color',colour,...
            'LineWidth',lineWidth);
    end
else
    error('Unknown plot style.');
end