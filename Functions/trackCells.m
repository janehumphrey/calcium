%%
% NAME: TRACK CELLS
% AUTHOR: JANE HUMPHREY (janehumphrey@outlook.com)

function [coord,trackArray,nCells,logIndex,cellIndex,firstFrame,finalFrame,landingFrame,xPos,yPos,cumDist,speed,...
    Results] = trackCells(stack,minInt,radiusPix,nFrames,Info,memory)

if nargin<5
    error('Not enough input arguments.');
end

coord = cell(nFrames,1);
peaks = cell(nFrames,1);
for iFrame = 1:nFrames
    coord{iFrame} = pkfnd(stack(:,:,iFrame),minInt,radiusPix*1.5);
    nPeaks = size(coord{iFrame},1);
    peaks{iFrame} = NaN(nPeaks,3);
    peaks{iFrame}(:,1:2) = coord{iFrame};
    peaks{iFrame}(:,3) = iFrame;
end
peaks = vertcat(peaks{:});

if nargin<6||isempty(memory)
    trackingParam.mem = 5;
else
    trackingParam.mem = memory;
end
trackingParam.dim = 2;
minLengthFrames = Info.minLength/Info.interval;
trackingParam.good = minLengthFrames*0.8; 
trackingParam.quiet = 1;
trackMat = track(peaks,radiusPix,trackingParam);

nTracks = trackMat(end,4);
goodTracks = true(nTracks,1);
trackArray = cell(nTracks,1);
firstFrame = NaN(nTracks,1);
finalFrame = NaN(nTracks,1);
lengthFrames = NaN(nTracks,1);
for iTrack = 1:nTracks
    trackArray{iTrack} = trackMat(trackMat(:,4)==iTrack,:);
    firstFrame(iTrack) = trackArray{iTrack}(1,3);
    finalFrame(iTrack) = trackArray{iTrack}(end,3);
    lengthFrames(iTrack) = finalFrame(iTrack)-firstFrame(iTrack);
    if lengthFrames(iTrack)<=minLengthFrames
        goodTracks(iTrack) = false;
    end
    if Info.earlyCells==0&&firstFrame(iTrack)<=10
        goodTracks(iTrack) = false;
    end
end

nCells = sum(goodTracks);
if nCells==0
    error('No cells selected.');
end
trackArray(~goodTracks) = [];
firstFrame(~goodTracks) = [];
finalFrame(~goodTracks) = [];
lengthFrames(~goodTracks) = [];
Results.length = lengthFrames*Info.interval;

logIndex = false(nFrames,nCells);
fullIndex = false(nFrames,nCells);
shortIndex = false(nFrames,nCells);
cellIndex = true(nCells,1);
xPos = NaN(nFrames,nCells);
yPos = NaN(nFrames,nCells);
landingFrame = NaN(nCells,1);
Results.xFinal = NaN(nCells,1);
Results.yFinal = NaN(nCells,1);
Results.attachment = NaN(nCells,1);
Results.landing = NaN(nCells,1);
Results.distance = NaN(nCells,1);
Results.speed = NaN(nCells,1);

for iCell = 1:nCells
    index = trackArray{iCell}(:,3);
    logIndex(index,iCell) = true;
    fullIndex(firstFrame(iCell):finalFrame(iCell),iCell) = true;
    shortIndex(firstFrame(iCell)+1:finalFrame(iCell)-1,iCell) = true;
    xPos(index,iCell) = trackArray{iCell}(:,1);
    yPos(index,iCell) = trackArray{iCell}(:,2);
    
    xLastValues = trackArray{iCell}(end-9:end,1);
    yLastValues = trackArray{iCell}(end-9:end,2);
    Results.xFinal(iCell) = round(mean(xLastValues));     
    Results.yFinal(iCell) = round(mean(yLastValues));
end

xPosSmooth = smoothdata(xPos,1,'gaussian',Info.smoothing*10/Info.interval);
yPosSmooth = smoothdata(yPos,1,'gaussian',Info.smoothing*10/Info.interval);
xPosSmooth(~fullIndex) = NaN;
yPosSmooth(~fullIndex) = NaN;
xDist = zeros(nFrames,nCells);
xDist(2:end,:) = diff(xPosSmooth,1,1)*Info.pixelSize;
yDist = zeros(nFrames,nCells);
yDist(2:end,:) = diff(yPosSmooth,1,1)*Info.pixelSize;
totalDist = (xDist.^2+yDist.^2).^0.5;
cumDist = cumsum(totalDist,'omitnan');
cumDist(~fullIndex) = NaN;
speed = gradient(cumDist')';
speed(~shortIndex) = NaN;
stopped = double(speed<Info.maxSpeed);
stopped(~shortIndex) = NaN;
changes = NaN(nFrames,nCells);
changes(2:end,:) = diff(stopped);

for iCell = 1:nCells
    cellStopped = stopped(shortIndex(:,iCell),iCell);
    if cellStopped(end)==true
        Results.attachment(iCell) = 1;
        if min(cellStopped)==1
            landingFrame(iCell) = firstFrame(iCell);
            Results.landing(iCell) = 0;
        else
            landingFrame(iCell) = find(changes(:,iCell)==1,1,'last');
            if landingFrame(iCell)>finalFrame(iCell)-50/Info.interval
                Results.attachment(iCell) = 0;
                landingFrame(iCell) = NaN;
            end
            Results.landing(iCell) = (landingFrame(iCell)-firstFrame(iCell))*Info.interval;
        end
    else
        Results.attachment(iCell) = 0;
    end
    cellDist = cumDist(:,iCell);
    cellDist(isnan(cellDist)) = [];
    Results.distance(iCell) = cellDist(end);
    Results.speed(iCell) = Results.distance(iCell)/Results.length(iCell);
end