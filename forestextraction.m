dataFolder = fullfile(tempdir,"forestData",filesep);
dataFile = dataFolder + "forestData.laz";
folderExists = exist(dataFolder,'dir');
fileExists = exist(dataFile,'file');
if ~folderExists
    mkdir(dataFolder);
end
if ~fileExists
    unzip('forestData.zip',dataFolder);
end
lasReader = lasFileReader(dataFile);
[ptCloud,pointAttributes] = readPointCloud(lasReader,"Attributes","ScanAngle");
figure
pcshow(ptCloud.Location)
title("Input Point Cloud")
% Segment Ground and extract non-ground and ground points
groundPtsIdx = segmentGroundSMRF(ptCloud);
nonGroundPtCloud = select(ptCloud,~groundPtsIdx);
groundPtCloud = select(ptCloud,groundPtsIdx);
% Visualize non-ground and ground points in magenta and green, respectively
figure
pcshowpair(nonGroundPtCloud,groundPtCloud)
title("Segmented Non-Ground and Ground Points")
groundPoints = groundPtCloud.Location;
% Eliminate duplicate points along x- and y-axes
[uniqueZ,uniqueXY] = groupsummary(groundPoints(:,3),groundPoints(:,1:2),@mean);
uniqueXY = [uniqueXY{:}];
% Create interpolant and use it to estimate ground elevation at each point
F = scatteredInterpolant(double(uniqueXY),double(uniqueZ),"natural");
estElevation = F(double(ptCloud.Location(:,1)),double(ptCloud.Location(:,2)));
% Normalize elevation by ground
normalizedPoints = ptCloud.Location;
normalizedPoints(:,3) = normalizedPoints(:,3) - estElevation;
% Visualize normalized points
figure
pcshow(normalizedPoints)
title("Point Cloud with Normalized Elevation")
% Set grid size to 10 meters per pixel and cutOffHeight to 2 meters
gridSize = 10;
cutOffHeight = 2;
leafAngDistribution = 0.5;
% Extract forest metrics
[canopyCover,gapFraction,leafAreaIndex] = helperExtractForestMetrics(normalizedPoints, ...
    pointAttributes.ScanAngle,gridSize,cutOffHeight,leafAngDistribution);
% Visualize forest metrics
hForestMetrics = figure;
axCC = subplot(2,2,1,Parent=hForestMetrics);
axCC.Position = [0.05 0.51 0.4 0.4];
imagesc(canopyCover,Parent=axCC)
title(axCC,"Canopy Cover")
axis off
colormap(gray)
axGF = subplot(2,2,2,Parent=hForestMetrics);
axGF.Position = [0.55 0.51 0.4 0.4];
imagesc(gapFraction,'Parent',axGF)
title(axGF,"Gap Fraction")
axis off
colormap(gray)
axLAI = subplot(2,2,[3 4],Parent=hForestMetrics);
axLAI.Position = [0.3 0.02 0.4 0.4];
imagesc(leafAreaIndex,Parent=axLAI)
title(axLAI,"Leaf Area Index")
axis off
colormap(gray)

% Set grid size to 0.5 meters per pixel 
gridRes = 0.5;
% Generate CHM
canopyModel = pc2dem(pointCloud(normalizedPoints),gridRes,CornerFillMethod="max");
% Clip invalid and negative CHM values to zero
canopyModel(isnan(canopyModel) | canopyModel<0) = 0;
% Perform gaussian smoothing to remove noise effects
H = fspecial("gaussian",[5 5],1);
canopyModel = imfilter(canopyModel,H,'replicate','same');
% Visualize CHM
figure
imagesc(canopyModel)
title('Canopy Height Model')
axis off
colormap(gray)

% Set minTreeHeight to 5 m 
minTreeHeight = 5;
% Detect tree tops
[treeTopRowId,treeTopColId] = helperDetectTreeTops(canopyModel,gridRes,minTreeHeight);
% Visualize treetops
figure
imagesc(canopyModel)
hold on
plot(treeTopColId,treeTopRowId,"rx",MarkerSize=3)
title("CHM with Detected Tree Tops")
axis off
colormap("gray")

% Segment individual trees
label2D = helperSegmentTrees(canopyModel,treeTopRowId,treeTopColId,minTreeHeight);
% Identify row and column id of each point in label2D and transfer labels
% to each point
rowId = ceil((ptCloud.Location(:,2) - ptCloud.YLimits(1))/gridRes) + 1;
colId = ceil((ptCloud.Location(:,1) - ptCloud.XLimits(1))/gridRes) + 1;
ind = sub2ind(size(label2D),rowId,colId);
label3D = label2D(ind);
% Extract valid labels and corresponding points
validSegIds = label3D ~= 0;
ptVeg = select(ptCloud,validSegIds);
veglabel3D = label3D(validSegIds);
% Assign color to each label
numColors = max(veglabel3D);
colorMap = randi([0 255],numColors,3)/255;
labelColors = label2rgb(veglabel3D,colorMap,OutputFormat="triplets");
% Visualize tree segments
figure
pcshow(ptVeg.Location,labelColors)
title("Individual Tree Segments")
view(2)

% Extract tree attributes
treeMetrics = helperExtractTreeMetrics(normalizedPoints,label3D);
% Display first 5 tree segments metrics
disp(head(treeMetrics,5));



function [canopyCover, gapFraction, leafAreaIndex] = helperExtractForestMetrics(normalizedPoints, ...
    scanAngles, gridSize, cutoffHeight,leafAngDistribution)

% Copyright 2021 The MathWorks, Inc.

% Marks points with values less than cutoffHeight as ground
groundPtsIdx = normalizedPoints(:,3) < cutoffHeight;

% Identify 2D index for each 3D point
xmin = min(normalizedPoints(:,1));
xmax = max(normalizedPoints(:,1));
ymin = min(normalizedPoints(:,2));
ymax = max(normalizedPoints(:,2));
nrows = round((ymax - ymin)/gridSize) + 1;
ncols = round((xmax - xmin)/gridSize) + 1;
rowId = round((normalizedPoints(:,2) - ymin)/gridSize) + 1;
colId = round((normalizedPoints(:,1) - xmin)/gridSize) + 1;
ind = sub2ind([nrows, ncols], rowId, colId);
[sortedInd, ptInd] = sort(ind);

% Identify the number of points at each cell
idCounts = histc(ind, unique(ind)); %#ok<HISTC>

% Initialize canpoy cover, gap fraction and leaf area index
canopyCover = zeros(nrows, ncols);
gapFraction = zeros(nrows, ncols);
leafAreaIndex = zeros(nrows, ncols);
endIdx = 0;

% Loop over valid cells
for i = 1: length(idCounts)
    gridIdx = endIdx + 1;
    endIdx = endIdx + idCounts(i);
    classes = groundPtsIdx(ptInd(gridIdx:endIdx));
    vegPts = sum(~classes);
    % CC = numVegPoints/totalNumPoints
    canopyCover(sortedInd(gridIdx)) = vegPts/idCounts(i);
    % GF = numGroundPoints/totalNumPoints
    gapFraction(sortedInd(gridIdx)) = 1 - canopyCover(sortedInd(gridIdx));
    %LAI = -0.5*cos(meanScanAngle)*log(GF)
    meanAng = sum(scanAngles(ptInd(gridIdx:endIdx)))/idCounts(i);
    leafAreaIndex(sortedInd(gridIdx)) = -(cosd(meanAng)*log(gapFraction(sortedInd(gridIdx))))/leafAngDistribution;
end
end

function [treeTopRowId,treeTopColId] = helperDetectTreeTops(canopyModel,gridRes,minTreeHeight)

% Copyright 2021 The MathWorks, Inc.

% Compute crown Radius and variable window Radius
crownRadius = (1.2 + 0.16 * canopyModel)/2;
windowRadius = max(round(crownRadius/gridRes),1);

% Mark window radius as 0 for points with elevation less than minTreeHeight 
windowRadius(canopyModel < minTreeHeight) = 0;
uniqueWindowRadius = sort(unique(windowRadius),'descend');

% Initialize non-canopy points to true and tree tops to false
nonCanopyPoints = true(size(canopyModel));
treeTopPoints = false(size(canopyModel));

% Loop over each unique radius and detect tree tops
for i=1:length(uniqueWindowRadius)-1
    % Create structuring element
    r = double(uniqueWindowRadius(i));
    SE = bwdist(padarray(1,[r r]))<=r;
    SE(ceil(size(SE,1)/2),ceil(size(SE,2)/2)) = 0;
    % Identify tree top ids
    treeTopIds = canopyModel>=imdilate(canopyModel,SE) & windowRadius==r & nonCanopyPoints;
    nonCanopyPoints = nonCanopyPoints & ~(bwdist(treeTopIds)<=r);
    treeTopPoints = treeTopPoints|treeTopIds;
end

% Identify row and column ids of tree top points
[treeTopRowId,treeTopColId] = find(treeTopPoints);
end

function label2D = helperSegmentTrees(canopyModel,treeTopRowId,treeTopColId,minTreeHeight)

% Copyright 2021 The MathWorks, Inc.

% Generate marker image 
markerImage = false(size(canopyModel));
vaildTreeTops = sub2ind(size(canopyModel),treeTopRowId,treeTopColId);
markerImage(vaildTreeTops) = true;

% Filter complement of CHM by minima imposition
canopyComplement = -canopyModel;
minImage = imimposemin(canopyComplement, markerImage);
label2D = watershed(minImage, 8);

% Remove labels for points with value less than minTreeHeight
label2D(canopyModel < minTreeHeight) = 0;

% Regroup labels
label2D = bwlabel(label2D ~= 0, 8);
end