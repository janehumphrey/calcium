%%
% NAME: CALCIUM
% AUTHOR: JANE HUMPHREY (janehumphrey@outlook.com)

tic;
close all;
clear variables;
clc;

% USER INPUT.
inputDir = '';  % Folder containing a single .tif file.
bgFile = '';  % Background file (optional).
resultsDir = '';  % Results folder.
Info.cellDiam = 10;  % Cell diameter, in um.
Info.pixelSize = 16/20;  % Pixel size, in um.
Info.interval = 1;  % Time between frames, in s.
Info.minLength = 200;  % Minimum track length, in s.
Info.minInt = 0;  % Intensity threshold for cell selection (SDs above mean). After running check Peaks.tif to adjust.
Info.smoothing = 5;  % Standard deviation for smoothing of traces, in s.
Info.maxSpeed = 0.1;  % Speed beneath which a cell is classed as stationary, in um/s.
Info.minGrad = 0.025;  % Threshold (intensity gradient) defining a calcium spike, in /s. After running check intensity
% traces to adjust.
Info.minHeight = 1.5;  % Minimum intensity of a calcium spike, relative to baseline (optional).
Info.minWidth = 0;  % Minimum duration of a calcium spike, in s (optional).
Info.hotHeight = 1.5;  % Intensity (at start of track) above which cells are classed as "coming in hot".
Info.earlyCells = true;  % If set to false, cells visible at the start are excluded.
%

functionDir = [pwd,filesep,'Functions'];
if isfolder(functionDir)
    addpath(genpath(functionDir));
    addpath(genpath([pwd,filesep,'Extra_functions']));
else
    addpath(genpath([fileparts(pwd),filesep,'Functions']));
    addpath(genpath([fileparts(pwd),filesep,'Extra_functions']));
end
if ~isfolder(inputDir)
    error('Directory doesn''t exist.');
end

[fileList,fileExt,nFiles] = selectFiles(inputDir,'tif');
if nFiles==1
    file = [inputDir,filesep,fileList{1},fileExt{1}];
else
    error('Too many files in directory.');
end

resultsDir = createDir(resultsDir,inputDir,'Calcium');
trueDir = [resultsDir,filesep,'Calcium_release'];
mkdir(trueDir);
falseDir = [resultsDir,filesep,'No_calcium_release'];
mkdir(falseDir);
hotDir = [resultsDir,filesep,'Came_in_hot'];
mkdir(hotDir);
movDir = [resultsDir,filesep,'Movement'];
mkdir(movDir);

radiusPix = round(Info.cellDiam/Info.pixelSize/2);
[fovRaw,width,height,nFrames] = readStack(file);
fovDp = double(fovRaw);
if isempty(bgFile)
    bgMean = min(fovDp(:));
else
    bgRaw = readStack(bgFile);
    bgMean = mean(bgRaw,3);
end
fovNoBg = fovDp-bgMean;

fovAvg = mean(fovNoBg,3);
fovFilt = imgaussfilt(fovAvg,radiusPix*10);
missing2d = fovFilt<=1;
missing3d = repmat(missing2d,[1,1,nFrames]);
fovFlat = fovNoBg./fovFilt;
fovFlat(missing3d) = 0;
gaussSd = round(radiusPix/2);
fovSmooth = imgaussfilt(fovFlat,gaussSd);

fovMean = mean(fovSmooth(:));
fovSd = mean(fovSmooth(:));
peakThresh = fovMean+Info.minInt*fovSd;
frameNos = (1:nFrames)';
time = frameNos*Info.interval;

[coord,trackArray,nCells,logIndex,cellIndex,firstFrame,finalFrame,landingFrame,xPos,yPos,cumDist,speed,Results] = ...
    trackCells(fovSmooth,peakThresh,radiusPix,nFrames,Info);

cellNos = (1:nCells)';
selectedCells = cellNos(cellIndex);

normTraces = NaN(nFrames,nCells);
gradTraces = NaN(nFrames,nCells);
spikeLoc = cell(nCells,1);
areas = cell(nCells,1);
cellId = cell(nCells,1);
Spikes.no = cell(nCells,1);
Spikes.time = cell(nCells,1);
Spikes.lag = cell(nCells,1);
Spikes.height = cell(nCells,1);
Spikes.width = cell(nCells,1);
Spikes.area = cell(nCells,1);
Results.caRelease = NaN(nCells,1);
Results.hot = NaN(nCells,1);
Results.nSpikes = NaN(nCells,1);
Results.initialTime = NaN(nCells,1);
Results.initialLag = NaN(nCells,1);
Results.maxHeight = NaN(nCells,1);
Results.maxWidth = NaN(nCells,1);
Results.combArea = NaN(nCells,1);

se = strel('disk',radiusPix,0);
warning('off','signal:findpeaks:largeMinPeakHeight');

for iCell = selectedCells'
    goodSpikes = [];
    index = trackArray{iCell}(:,3);
    intTrace = NaN(nFrames,1);
    for iFrame = index'
        roiX = xPos(iFrame,iCell)-radiusPix:xPos(iFrame,iCell)+radiusPix;
        roiY = yPos(iFrame,iCell)-radiusPix:yPos(iFrame,iCell)+radiusPix;
        goodWidth = roiX>0&roiX<=width;
        goodHeight = roiY>0&roiY<=height;
        roiBox = fovSmooth(roiY(goodHeight),roiX(goodWidth),iFrame);
        roiDisk = roiBox.*se.Neighborhood(goodHeight,goodWidth);
        intTrace(iFrame) = nanmean(roiDisk(:));
    end
    
    smoothTrace = smoothdata(intTrace,'gaussian',Info.smoothing*5/Info.interval);
    smoothTrace(1:firstFrame(iCell)-1) = NaN;
    smoothTrace(finalFrame(iCell)+1:end) = NaN;
    tempSmooth = smoothTrace;
    tempSmooth(isnan(tempSmooth)) = [];
    sortedTrace = sort(tempSmooth);
    cutoff = round(size(sortedTrace,1)/4);
    refSmooth = sortedTrace(1:cutoff);
    baselineEst = mean(refSmooth);
    normTraces(:,iCell) = smoothTrace/baselineEst;
    gradTraces(:,iCell) = gradient(normTraces(:,iCell));
    tempGrad = gradTraces(:,iCell);
    [peakHeight,peakLoc,~,peakProm] = findpeaks(tempGrad,'minPeakHeight',Info.minGrad);
    spikeStart = peakLoc(peakProm>0.75*peakHeight);
    spikeStart(spikeStart<firstFrame(iCell)+20/Info.interval) = [];
    if ~isempty(spikeStart)
        heightAtStart = normTraces(spikeStart,iCell);
        nSpikes = size(spikeStart,1);
        overlapping = false(nSpikes,1);
        spikeEnd = NaN(nSpikes,1);
        spikeInt = NaN(nFrames,nSpikes);
        for iSpike = 1:nSpikes
            if max(spikeStart(iSpike)<spikeEnd(:))
                overlapping(iSpike) = true;
            end
            afterSpikeStart = frameNos>spikeStart(iSpike);
            belowStartHeight = normTraces(:,iCell)<heightAtStart(iSpike);
            if ~isempty(find(afterSpikeStart&belowStartHeight,1))
                spikeEnd(iSpike) = find(afterSpikeStart&belowStartHeight,1);
            else
                spikeEnd(iSpike) = nFrames;
            end
            spikeInt(spikeStart(iSpike):spikeEnd(iSpike),iSpike) = normTraces(spikeStart(iSpike):spikeEnd(iSpike),...
                iCell);
        end
        spikeHeight = max(spikeInt,[],1)';
        spikeWidth = (spikeEnd-spikeStart)*(Info.interval);
        goodHeight = spikeHeight>=Info.minHeight;
        goodWidth = spikeWidth>=Info.minWidth;
        goodSpikes = goodHeight&goodWidth&~overlapping;
    end
    heightComingIn = max(normTraces(firstFrame(iCell):(firstFrame(iCell)+round(25/Info.interval)),iCell));
    if any(goodSpikes)==1
        spikeLoc{iCell} = spikeStart(goodSpikes);
        Spikes.time{iCell} = (spikeStart(goodSpikes)-firstFrame(iCell))*Info.interval;
        Spikes.lag{iCell} = (spikeStart(goodSpikes)-landingFrame(iCell))*Info.interval;
        if isnan(landingFrame(iCell))
            Spikes.lag{iCell}(:) = -nFrames;
        end
        Spikes.height{iCell} = spikeHeight(goodSpikes);
        Spikes.width{iCell} = spikeWidth(goodSpikes);
        areas{iCell} = spikeInt(:,goodSpikes);
        Spikes.area{iCell} = nansum(areas{iCell}-1,1)';
        
        Results.caRelease(iCell) = 1;
        Results.hot(iCell) = 0;
        Results.nSpikes(iCell) = size(spikeLoc{iCell},1);
        Results.initialTime(iCell) = Spikes.time{iCell}(1);
        Results.initialLag(iCell) = Spikes.lag{iCell}(1);
        Results.maxHeight(iCell) = max(Spikes.height{iCell});
        Results.maxWidth(iCell) = max(Spikes.width{iCell});
        Results.combArea(iCell) = sum(Spikes.area{iCell});
        
        cellId{iCell} = ones(Results.nSpikes(iCell),1)*iCell;
        Spikes.no{iCell} = (1:Results.nSpikes(iCell))';
    elseif heightComingIn>Info.hotHeight
        Results.hot(iCell) = 1;
    else
        Results.caRelease(iCell) = 0;
        Results.hot(iCell) = 0;
    end
end

cellsTrue = cellNos(Results.caRelease==1);
cellsFalse = cellNos(Results.caRelease==0);
cellsHot = cellNos(Results.hot==1);
maxDigits = size(num2str(cellNos(end)),2);
format = ['%0.',num2str(maxDigits),'u'];
cellList = num2str(cellNos,format);

prepareFigure(width/height,width/height);
ax = prepareAxes('image',[],[],[1-10,width+10],[1-10,height+10]);
xCorners = [1-10,width+10,width+10,1-10]';
yCorners = [1-10,1-10,height+10,height+10]';
lastFrame = fovSmooth(:,:,end);
minInt = quantile(lastFrame(:),0.2);    
maxInt = quantile(lastFrame(:),0.99);    
lutLimits = [minInt,maxInt];
if radiusPix/height>0.01
    markerSize = 100;
    fontSize = 11;
else
    markerSize = 60;
    fontSize = 8;
end
peaksFile = [resultsDir,filesep,'Peaks.tif'];
pointsFile = [resultsDir,filesep,'Cells.tif'];
labelsFile = [resultsDir,filesep,'Labelled_cells.tif'];

for iFrame = 100:100:nFrames
    patch(xCorners,yCorners,[0,0.5,1],'EdgeColor','none');
    imagesc(fovSmooth(:,:,iFrame),lutLimits);
    plotData('points',coord{iFrame}(:,1),coord{iFrame}(:,2),[0,0.5,1],[],[],[],'none',markerSize);
    export_fig(peaksFile,'-tif','-append');
    ax.Children.delete;
end

for iFrame = 10:10:nFrames
    patch(xCorners,yCorners,[1,0.5,0],'EdgeColor','none');
    imagesc(fovSmooth(:,:,iFrame),lutLimits);

    trackIndex = logIndex(iFrame,:)';
    trackList = cellNos(trackIndex);
    tracksTrue = intersect(cellsTrue,trackList);
    tracksFalse = intersect(cellsFalse,trackList);
    tracksHot = intersect(cellsHot,trackList);
    xTrue = xPos(iFrame,tracksTrue)';
    yTrue = yPos(iFrame,tracksTrue)';
    xFalse = xPos(iFrame,tracksFalse)';
    yFalse = yPos(iFrame,tracksFalse)';
    xHot = xPos(iFrame,tracksHot)';
    yHot = yPos(iFrame,tracksHot)';
    trueStr = num2str(tracksTrue);
    falseStr = num2str(tracksFalse);
    hotStr = num2str(tracksHot);

    trueArray = cellstr(trueStr);
    falseArray = cellstr(falseStr);
    hotArray = cellstr(hotStr);
    
    plotData('points',xTrue,yTrue,[],[],[],[],'none',markerSize);
    plotData('points',xFalse,yFalse,[1,0,0],[],[],[],'none',markerSize);
    plotData('points',xHot,yHot,[1,1,0],[],[],[],'none',markerSize);
    export_fig(pointsFile,'-tif','-append');
    ax.Children(1:end-2).delete;
    
    plotData('text',xTrue,yTrue,[],[],[],[],[],[],trueArray,fontSize);
    plotData('text',xFalse,yFalse,[1,0,0],[],[],[],[],[],falseArray,fontSize);
    plotData('text',xHot,yHot,[1,1,0],[],[],[],[],[],hotArray,fontSize);
    export_fig(labelsFile,'-tif','-append');
    ax.Children.delete;
end
close;

prepareFigure(1);
prepareAxes('image',[],[],[1,width],[1,height]);
for iCell = 1:nCells
    if Results.caRelease(iCell)==1
        colour = [1,0.5,0];
    elseif Results.caRelease(iCell)==0
        colour = [1,0,0];
    elseif Results.hot(iCell)==1
        colour = [1,1,0];
    else
        colour = [0,0.5,1];
    end
    intIndex = logIndex(:,iCell);
    plotData('line',xPos(intIndex,iCell),yPos(intIndex,iCell),colour);
end
tracksBase = [resultsDir,filesep,'Tracks'];
saveImage(tracksBase);
close;

prepareFigure(1);
prepareAxes('image',[],[],[1,width],[1,height]);
colormap('jet');
scrambledCells = randperm(nCells);
for iCell = scrambledCells
    intIndex = logIndex(:,iCell);
    xValues = xPos(intIndex,iCell);
    yValues = yPos(intIndex,iCell);
    zValues = zeros(sum(intIndex),1);
    colorValues = time(intIndex);
    surface([xValues,xValues],[yValues,yValues],[zValues,zValues],[colorValues,colorValues],'facecol','no','edgecol',...
        'interp','linew',3);
end
legendHandle = colorbar('Location','southoutside','AxisLocation','in','Ticks',[],'FontSize',15,...
    'TickLabelInterpreter','tex');
legendHandle.Label.String = 'time \rightarrow';
tracksBase = [resultsDir,filesep,'Tracks_wrt_time'];
saveImage(tracksBase);
close;

nHot = size(cellsHot,1);
binEdges = [-Inf;time];
includedCells = (nCells-nHot);
ordered = sort(Results.initialTime,1);
binned = cumsum(histcounts(ordered,binEdges))';
triggerTime = binned/includedCells*100;
ordered = sort(Results.initialLag,1);
binned = cumsum(histcounts(ordered,binEdges))';
triggerLag = binned/includedCells*100;

prepareFigure;
prepareAxes('plot','time (s)','triggered cells (%)',[0,nFrames],[0,100]);
[~,index] = unique(triggerTime,'stable');
timeX = time([index;nFrames]);
timeY = triggerTime([index;nFrames]);
plotData('line',timeX,timeY);
timeBase = [resultsDir,filesep,'Triggering'];
saveImage(timeBase);
close;

prepareFigure;
prepareAxes('plot','time wrt landing (s)','triggered cells (%)',[0,nFrames],[0,100]);
[~,index] = unique(triggerLag,'stable');
lagX = time([index;nFrames]);
lagY = triggerLag([index;nFrames]);
plotData('line',lagX,lagY);
lagBase = [resultsDir,filesep,'Triggering_wrt_landing'];
saveImage(lagBase);
close;

if ~isempty(selectedCells)
    prepareFigure(0.75);
    subplot(2,1,1);
    normMax = max(normTraces(:));
    if normMax>nanmean(normTraces(:))*10
        normMax = nanmean(normTraces(:))*10;
    end
    topAx = prepareAxes('plot','time (s)','normalised intensity',[0,time(end)],[0,normMax]);
    subplot(2,1,2);
    gradMax = max(gradTraces(:));
    if gradMax>Info.minGrad*10
        gradMax = Info.minGrad*10;
    elseif gradMax<Info.minGrad*2
        gradMax = Info.minGrad*2;
    end
    gradLabel = ['intensity gradient (',getUnit('perS'),')'];
    bottomAx = prepareAxes('plot','time (s)',gradLabel,[0,time(end)],[-gradMax,gradMax]);
    gradThresh = [Info.minGrad,Info.minGrad];
    plotData('line',[0,time(end)],gradThresh,[0,0,0],':',2);

    for iCell = selectedCells'
        subplot(2,1,1);
        if Results.caRelease(iCell)==1
            for iSpike = 1:Results.nSpikes(iCell)
                patches = area(time,areas{iCell}(:,iSpike),1,'LineStyle',':','LineWidth',2,'FaceColor',...
                    [0.85,0.85,0.85]);
                patches.BaseLine.Visible = 'off';
            end
        end
        if Results.caRelease(iCell)==1
            plotData('line',time,normTraces(:,iCell));
        elseif Results.caRelease(iCell)==0
            plotData('line',time,normTraces(:,iCell),[1,0,0]);
        elseif Results.hot(iCell)==1
            plotData('line',time,normTraces(:,iCell),[1,1,0]);
        end
        if ~isnan(landingFrame(iCell))
            landingTime = time(landingFrame(iCell));
            landingNorm = normTraces(landingFrame(iCell),iCell);
            plotData('points',landingTime,landingNorm,[0,0,0],[],[],'d',[],50);
        end
        subplot(2,1,2);
        if Results.caRelease(iCell)==1
            plotData('line',time,gradTraces(:,iCell));
        elseif Results.caRelease(iCell)==0
            plotData('line',time,gradTraces(:,iCell),[1,0,0]);
        elseif Results.hot(iCell)==1
            plotData('line',time,gradTraces(:,iCell),[1,1,0]);
        end
        if Results.caRelease(iCell)==1
            spikeTime = time(spikeLoc{iCell});
            spikeGrad = gradTraces(spikeLoc{iCell},iCell);
            plotData('points',spikeTime,spikeGrad,[0,0,0],[],[],[],[],50);
            imageBase = [trueDir,filesep,'Cell_',cellList(iCell,:)];
        elseif Results.caRelease(iCell)==0
            imageBase = [falseDir,filesep,'Cell_',cellList(iCell,:)];
        elseif Results.hot(iCell)==1
            imageBase = [hotDir,filesep,'Cell_',cellList(iCell,:)];
        end
        saveImage(imageBase);
        topAx.Children.delete;
        bottomAx.Children(1:end-1).delete;
    end
    close;
end

[landingX,landingY] = plotMovement(time,nFrames,nCells,cumDist,speed,landingFrame,resultsDir,movDir,cellList,Info.maxSpeed,Results);

excel = actxserver('Excel.Application');

intBase = repmat('Cell ',nCells,1);
intMat = [intBase,num2str(cellNos)];
intTitles = cellstr(intMat)';
intValues = num2cell(normTraces);
intText = [{'Time (s)'},intTitles;num2cell(time),intValues];
intFile = [resultsDir,filesep,'Traces.xlsx'];
createSpreadsheet(excel,intFile,intText);

timeText = [{'Time','Triggered cells (%)'};num2cell(timeX),num2cell(timeY)];
timeFile = [resultsDir,filesep,'Triggering.xlsx'];
createSpreadsheet(excel,timeFile,timeText,0);

lagText = [{'Time','Triggered cells (%)'};num2cell(lagX),num2cell(lagY)];
lagFile = [resultsDir,filesep,'Triggering_wrt_landing.xlsx'];
createSpreadsheet(excel,lagFile,lagText,0);
    
landingText = [{'Time','Attached cells (%)'};num2cell(landingX),num2cell(landingY)];
landingFile = [resultsDir,filesep,'Landing.xlsx'];
createSpreadsheet(excel,landingFile,landingText,0);

cellDiamStr = ['Cell diameter (',getUnit('um'),')'];
pixelSizeStr = ['Pixel size (',getUnit('um'),')'];
maxSpeedStr = ['Max speed (',getUnit('umPerS'),')'];
minGradStr = ['Min gradient (',getUnit('perS'),')'];
infoValues = struct2cell(Info)';
infoText = [{cellDiamStr,pixelSizeStr,'Interval (s)','Min track length (s)','Intensity threshold (SDs)',...
    'Smoothing (SDs)',maxSpeedStr,minGradStr,'Min spike height','Min spike duration (s)',...
    'Min height for "coming in hot"','Early cells included'};infoValues];
infoFile = [resultsDir,filesep,'Info.xlsx'];
createSpreadsheet(excel,infoFile,infoText);

cellId = vertcat(cellId{:});
spikesFields = fieldnames(Spikes);
for iField = 1:size(spikesFields,1)
    fieldName = spikesFields{iField};
    Spikes.(fieldName) = vertcat(Spikes.(fieldName){:});
    Spikes.(fieldName) = round(Spikes.(fieldName),3,'significant');
end
spikesValues = cell2mat(struct2cell(Spikes)');
spikesMean = round(nanmean(spikesValues,1),3,'significant');
spikesSd = round(nanstd(spikesValues,0,1),3,'significant');
if ~isempty(spikesValues)
    spikesHeadings = {'Cell','Spike','Time wrt start of track (s)','Time wrt landing (s)','Spike height',...
        'Spike duration (s)','Integrated intensity'};
    spikesText = [spikesHeadings;num2cell(cellId),num2cell(spikesValues);{'Mean'},num2cell(spikesMean);{'SD'},...
        num2cell(spikesSd)];
    spikesFile = [resultsDir,filesep,'Spikes.xlsx'];
    createSpreadsheet(excel,spikesFile,spikesText,2);
end

distStr = ['Distance travelled (',getUnit('um'),')'];
speedStr = ['Speed (',getUnit('umPerS'),')'];
resultsFields = fieldnames(Results);
for iField = 1:size(resultsFields,1)
    fieldName = resultsFields{iField};
    Results.(fieldName) = round(Results.(fieldName),3,'significant');
end
resultsValues = cell2mat(struct2cell(Results)');
resultsMean = round(nanmean(resultsValues,1),3,'significant');
resultsSd = round(nanstd(resultsValues,0,1),3,'significant');
resultsHeadings = {'Cell','Track length (s)','Final X position','Final Y position','Attachment',...
    'Time of landing (s)',distStr,speedStr,'Calcium release','Came in hot','No of spikes',...
    'Time of first spike wrt start of track (s)','Time of first spike wrt landing (s)','Max spike height',...
    'Max spike duration (s)','Total integrated intensity (all spikes)'};
resultsText = [resultsHeadings;num2cell(cellNos),num2cell(resultsValues);{'Mean'},num2cell(resultsMean);{'SD'},...
    num2cell(resultsSd)];
resultsFile = [resultsDir,filesep,'Results.xlsx'];
createSpreadsheet(excel,resultsFile,resultsText,2);

excel.Quit;
excel.delete;

toc;