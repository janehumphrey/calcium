%%
% NAME: PLOT MOVEMENT
% AUTHOR: JANE HUMPHREY (janehumphrey@outlook.com)

function [landingX,landingY] = plotMovement(time,nFrames,nCells,cumDist,speed,landingFrame,resultsDir,movDir,...
    cellList,maxSpeed,Results)

if nargin<8
    error('Not enough input arguments.');
end

binEdges = [-Inf;time];
ordered = sort(Results.landing,1);
binned = cumsum(histcounts(ordered,binEdges))';
attachment = binned/nCells*100;

prepareFigure;
prepareAxes('plot','time (s)','attached cells (%)',[0,nFrames],[0,100]);
[~,index] = unique(attachment,'stable');
landingX = time([index;nFrames]);
landingY = attachment([index;nFrames]);
plotData('line',landingX,landingY);
landingBase = [resultsDir,filesep,'Landing'];
saveImage(landingBase);
close;

distMax = max(cumDist(:));
speedMax = max(speed(:));
distLabel = ['distance (',getUnit('um'),')'];
speedLabel = ['speed (',getUnit('umPerS'),')'];
prepareFigure(0.75);
subplot(2,1,1);
topAx = prepareAxes('plot','time (s)',distLabel,[0,time(end)],[0,distMax]);
subplot(2,1,2);
bottomAx = prepareAxes('plot','time (s)',speedLabel,[0,time(end)],[0,speedMax]);
speedThresh = [maxSpeed,maxSpeed];
plotData('line',[0,time(end)],speedThresh,[0,0,0],':',2);

for iCell = 1:nCells
    subplot(2,1,1);
    if ~isfield(Results,'caRelease')||Results.caRelease(iCell)==1
        colour = [1,0.5,0];
    elseif Results.caRelease(iCell)==0
        colour = [1,0,0];
    elseif Results.hot(iCell)==1
        colour = [1,1,0];
    else
        colour = [0,0.5,1];     
    end
    plotData('line',time,cumDist(:,iCell),colour);
    if ~isnan(landingFrame(iCell))
        landingTime = time(landingFrame(iCell));
        landingDist = cumDist(landingFrame(iCell),iCell);
        plotData('points',landingTime,landingDist,[0,0,0],[],[],'d',[],50);
    end
    subplot(2,1,2);
    plotData('line',time,speed(:,iCell),colour);
    if ~isnan(landingFrame(iCell))
        landingSpeed = speed(landingFrame(iCell),iCell);
        plotData('points',landingTime,landingSpeed,[0,0,0],[],[],'d',[],50);
    end
    movBase = [movDir,filesep,'Cell_',cellList(iCell,:)];
    saveImage(movBase);
    topAx.Children.delete;
    bottomAx.Children(1:end-1).delete;
end
close;