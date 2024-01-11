datafile = 'aerialMap.tar.gz';
wayPointsfile = 'gTruthWayPoints.mat';

% Generate a lidar scan at each waypoint using the helper function
[pClouds,orgMap,gTruthWayPts] = helperCreateDataset(datafile,wayPointsfile);

% Create a figure window to visualize the ground truth map and waypoints
hFigGT = figure;
axGT = axes(Parent=hFigGT,Color='black');

% Visualize the ground truth waypoints
pcshow(gTruthWayPts,'red',MarkerSize=150,Parent=axGT)
hold on

% Visualize the original map covered by the robot
pcshow(orgMap,MarkerSize=10,Parent=axGT)
hold off

% Customize the axis labels
axis on
xlabel(axGT,'X (m)')
ylabel(axGT,'Y (m)')
zlabel(axGT,'Z (m)')
title(axGT,'Ground Truth Map And Robot Trajectory')

% Specify limits for the player
xlimits = [-90 90];
ylimits = [-90 90];
zlimits = [-625 -587];

% Create a pcplayer object to visualize the lidar scans
lidarPlayer = pcplayer(xlimits,ylimits,zlimits);

% Customize the pcplayer axis labels
xlabel(lidarPlayer.Axes,'X (m)')
ylabel(lidarPlayer.Axes,'Y (m)')
zlabel(lidarPlayer.Axes,'Z (m)')
title(lidarPlayer.Axes,'Lidar Scans')

% Loop over and visualize the data
for l = 1:length(pClouds)

    % Extract the lidar scan
    ptCloud = pClouds(l);

    % Update the lidar display
    view(lidarPlayer,ptCloud)
    pause(0.05)
end

skipFrames = 3;
gridStep = 1.5; % in meters
neighbors = 60;
matchThreshold = 0.1;
matchRatio = 0.97;
maxDistance = 1;
maxNumTrails = 3000;
inlierRatio = 0.1;
alignGridStep = 1.2;
loopClosureSearchRadius = 7.9;
nScansPerSubmap = 3;
subMapThresh = 15;
loopClosureThreshold = 0.6;
rmseThreshold = 0.6;

pGraph = poseGraph3D;

% Default serialized upper-right triangle of a 6-by-6 Information Matrix
infoMat = [1 0 0 0 0 0 1 0 0 0 0 1 0 0 0 1 0 0 1 0 1];

% Allocate memory to store submaps
subMaps = cell(floor(length(pClouds)/(skipFrames*nScansPerSubmap)),1);

% Allocate memory to store each submap pose
subMapPoses = zeros(round(length(pClouds)/(skipFrames*nScansPerSubmap)),3);

% Initialize variables to store accepted scans and their transforms for
% submap creation
pcAccepted = pointCloud.empty(0);
tformAccepted = rigidtform3d.empty(0);

% Initialize variable to store relative poses from the feature-based approach
% without pose graph optimization
fpfhTform = rigidtform3d.empty(0);

% Counter to track the number of scans added
count = 1;

% Set to 1 to visualize processed lidar scans during build process
viewPC = 0;

% Create a pcplayer object to view the lidar scans while
% processing them sequentially, if viewPC is enabled
if viewPC == 1
    pplayer = pcplayer(xlimits,ylimits,zlimits);

    % Customize player axis labels
    xlabel(pplayer.Axes,'X (m)')
    ylabel(pplayer.Axes,'Y (m)')
    zlabel(pplayer.Axes,'Z (m)')
    title(pplayer.Axes,'Processed Scans')
end

% Create a figure window to visualize the estimated trajectory
hFigTrajUpdate = figure;
axTrajUpdate = axes(Parent=hFigTrajUpdate,Color='black');
title(axTrajUpdate,'Sensor Pose Trajectory')

rng('default') % Set random seed to guarantee consistent results in pcmatchfeatures
for FrameIdx = 1:skipFrames:length(pClouds)
     % Downsample the current scan
    curScan = pcdownsample(pClouds(FrameIdx),gridAverage=gridStep);
    if viewPC == 1

        % Visualize down sampled point cloud
        view(pplayer,curScan)
    end
 % Extract FPFH features
    curFeature = extractFPFHFeatures(curScan,NumNeighbors=neighbors);

    if FrameIdx == 1

        % Update the acceptance scan and its tform
        pcAccepted(count,1) = curScan;
        tformAccepted(count,1) = rigidtform3d;

        % Update the initial pose to the first waypoint of ground truth for
        % comparison
        fpfhTform(count,1) = rigidtform3d([0 0 0],gTruthWayPts(1,:));
    else

        % Identify correspondences by matching current scan to previous scan
        indexPairs = pcmatchfeatures(curFeature,prevFeature,curScan,prevScan, ...
            MatchThreshold=matchThreshold,RejectRatio=matchRatio);
        matchedPrevPts = select(prevScan,indexPairs(:,2));
        matchedCurPts = select(curScan,indexPairs(:,1));

        % Estimate relative pose between current scan and previous scan
        % using correspondences
        tform1 = estgeotform3d(matchedCurPts.Location, ...
            matchedPrevPts.Location,'rigid',MaxDistance=maxDistance, ...
            MaxNumTrials=maxNumTrails);

        % Perform ICP registration to fine-tune relative pose
        tform = pcregistericp(curScan,prevScan,InitialTransform=tform1, ...
            InlierRatio=inlierRatio);

        relPose = [tform2trvec(tform.A) tform2quat(tform.A)];

        % Add relative pose to pose graph
        addRelativePose(pGraph,relPose);

         % Update counter and store accepted scans and their poses
        count = count + 1;
        pcAccepted(count,1) = curScan;
        accumPose = pGraph.nodes(height(pGraph.nodes));
        tformAccepted(count,1) = rigidtform3d((trvec2tform(accumPose(1:3)) * ...
            quat2tform(accumPose(4:7))));

        % Update estimated poses
        fpfhTform(count) = rigidtform3d(fpfhTform(count-1).A * tform.A);
    end

     currSubMapId = floor(count/nScansPerSubmap);
    if rem(count,nScansPerSubmap) == 0

        % Align an array of lidar scans to create a submap
        subMaps{currSubMapId} = pcalign(...
            pcAccepted((count - nScansPerSubmap + 1):count,1), ...
            tformAccepted((count - nScansPerSubmap + 1):count,1), ...
            alignGridStep);

        % Assign center scan pose as pose of submap
        subMapPoses(currSubMapId,:) = tformAccepted(count - ...
            floor(nScansPerSubmap/2),1).Translation;
    end

 if currSubMapId > subMapThresh
        mostRecentScanCenter = pGraph.nodes(height(pGraph.nodes));

        % Estimate possible loop closure candidates by matching current
        % scan with submaps
        [loopSubmapIds,~] = helperEstimateLoopCandidates(subMaps,curScan, ...
            subMapPoses,mostRecentScanCenter,currSubMapId,subMapThresh, ...
            loopClosureSearchRadius,loopClosureThreshold);

        if ~isempty(loopSubmapIds)
            rmseMin = inf;

            % Estimate the best match for the current scan from the matching submap ids
            for k = 1:length(loopSubmapIds)

                % Check every scan within the submap
                for kf = 1:nScansPerSubmap
                    probableLoopCandidate = ...
                        loopSubmapIds(k)*nScansPerSubmap - kf + 1;
                    [pose_Tform,~,rmse] = pcregistericp(curScan, ...
                        pcAccepted(probableLoopCandidate));

                    % Update the best loop closure candidate
                    if rmse < rmseMin
                        rmseMin = rmse;
                        matchingNode = probableLoopCandidate;
                        Pose = [tform2trvec(pose_Tform.A) ...
                            tform2quat(pose_Tform.A)];
                         end
                end
            end

            % Check if loop closure candidate is valid
            if rmseMin < rmseThreshold

                % Add relative pose of loop closure candidate to pose graph
                addRelativePose(pGraph,Pose,infoMat,matchingNode, ...
                    pGraph.NumNodes);
            end
        end
    end

    % Update previous point cloud and feature
    prevScan = curScan;
    prevFeature = curFeature;

    % Visualize the estimated trajectory from the accepted scan.
    show(pGraph,IDs='off',Parent=axTrajUpdate);
    drawnow
end

pGraph = optimizePoseGraph(pGraph,FirstNodePose=[gTruthWayPts(1,:) 1 0 0 0]);

% Get estimated trajectory from pose graph
pGraphTraj = pGraph.nodes;

% Get estimated trajectory from feature-based registration without pose
% graph optimization
fpfhEstimatedTraj = zeros(count,3);
for i = 1:count
    fpfhEstimatedTraj(i,:) = fpfhTform(i).Translation;
end

% Create a figure window to visualize the ground truth and estimated
% trajectories
hFigTraj = figure;
axTraj = axes(Parent=hFigTraj,Color='black');
plot3(fpfhEstimatedTraj(:,1),fpfhEstimatedTraj(:,2),fpfhEstimatedTraj(:,3), ...
    'r*',Parent=axTraj)
hold on
plot3(pGraphTraj(:,1),pGraphTraj(:,2),pGraphTraj(:,3),'b.',Parent=axTraj)
plot3(gTruthWayPts(:,1),gTruthWayPts(:,2),gTruthWayPts(:,3),'go',Parent=axTraj)
hold off
axis equal
view(axTraj,2)
xlabel(axTraj,'X (m)')
ylabel(axTraj,'Y (m)')
zlabel(axTraj,'Z (m)')
title(axTraj,'Trajectory Comparison')
legend(axTraj,'Estimated trajectory without pose graph optimization', ...
    'Estimated trajectory with pose graph optimization', ...
    'Ground Truth Trajectory','Location','bestoutside')

% Get the estimated trajectory from poses
estimatedTraj = pGraphTraj(:,1:3);

% Convert relative poses to rigid transformations
estimatedTforms = rigidtform3d.empty(0);
for idx=1:pGraph.NumNodes
    pose = pGraph.nodes(idx);
    rigidPose = rigidtform3d((trvec2tform(pose(1:3)) * quat2tform(pose(4:7))));
    estimatedTforms(idx,1) = rigidPose;
end

% Create global map from processed point clouds and their relative poses
globalMap = pcalign(pcAccepted,estimatedTforms,alignGridStep);

% Create a figure window to visualize the estimated map and trajectory
hFigTrajMap = figure;
axTrajMap = axes(Parent=hFigTrajMap,Color='black');
pcshow(estimatedTraj,'red',MarkerSize=150,Parent=axTrajMap)
hold on
pcshow(globalMap,MarkerSize=10,Parent=axTrajMap)
hold off

% Customize axis labels
axis on
xlabel(axTrajMap,'X (m)')
ylabel(axTrajMap,'Y (m)')
zlabel(axTrajMap,'Z (m)')
title(axTrajMap,'Estimated Robot Trajectory And Generated Map')

% Create a figure window to display both the ground truth map and estimated map
hFigMap = figure(Position=[0 0 700 400]);
axMap1 = subplot(1,2,1,Color='black',Parent=hFigMap);
axMap1.Position = [0.08 0.2 0.4 0.55];
pcshow(orgMap,Parent=axMap1)
axis on
xlabel(axMap1,'X (m)')
ylabel(axMap1,'Y (m)')
zlabel(axMap1,'Z (m)')
title(axMap1,'Ground Truth Map')

axMap2 = subplot(1,2,2,Color='black',Parent=hFigMap);
axMap2.Position = [0.56 0.2 0.4 0.55];
pcshow(globalMap,Parent=axMap2)
axis on
xlabel(axMap2,'X (m)')
ylabel(axMap2,'Y (m)')
zlabel(axMap2,'Z (m)')
title(axMap2,'Estimated Map')




function [pClouds, orgMap, gTruth] = helperCreateDataset(datafile, wayPointsfile)
% The function loads aerial data and ground truth waypoints from data file
% and create lidar scans at each waypoint.
    dataFolder = fullfile(tempdir, 'aerialMapData', filesep);
    aerialFile = dataFolder + "aerialMap.laz";
    
    % Check whether the folder and data file already exist or not
    folderExists = exist(dataFolder, 'dir');
    fileExists = exist(aerialFile, 'file');
    % Create a new folder if it is not exist
    if ~folderExists
        mkdir(dataFolder);
    end
    
    % Extract aerial data file if it is not exist
    if ~fileExists
        untar(datafile, dataFolder);
    end
  
    % The data is collected at an altitude of 700m
    altitude = 700;
    
    % Load laz file
    lazReader = lasFileReader(aerialFile);
    lazPtCloud = readPointCloud(lazReader);

    % Convert the location field to single datatype
    lazPtCloud = pointCloud(single(lazPtCloud.Location), ...
        Intensity=lazPtCloud.Intensity);
    
    % Transform the point cloud such that it will be centered around
    % origin. This helps to create fasteroccupancy map and easycomputation
    xshift = sum(lazPtCloud.XLimits)/2;
    yshift = sum(lazPtCloud.YLimits)/2;
    tform = rigidtform3d([0 0 0],[-xshift -yshift 0]);
    tformPt = pctransform(lazPtCloud, tform);

    % Create occupancy map for the data
    pose = [ 0 0 0 1 0 0 0];
    map3D = occupancyMap3D(1);
    insertPointCloud(map3D,pose,tformPt.Location,300);

    % Load Ground truth trajectory and adjust it w.r.t transformed data
    traj = load(wayPointsfile);
    gTruth = [traj.gTruthTraj(:,1)-xshift traj.gTruthTraj(:,2)-yshift repmat(altitude,height(traj.gTruthTraj),1)];
    
    % Create directional angles with horizontal FOV of -13 to 13 and
    % vertical FOV of -2.5 to 2.5
    [xx,yy]=meshgrid(-13:0.08:13, -2.5:0.2:2.5);
    dirAngles = [xx(:) yy(:) repmat(90,numel(xx),1)];

    % Extract lidar scans at each waypoint using rayIntesection
    mapPoints = [];
    pClouds = pointCloud.empty();
    for i = 1:height(gTruth)
        % calculate pose at each waypoint.
        if i == 1
            posex = -(0.5*pi);
        else
            posex = atan((gTruth(i,2)-gTruth(i-1,2))/(gTruth(i,1)-gTruth(i-1,1)));
        end
        rot = eul2quat([posex+0.5*pi 0 pi]);
        
        % Apply rayIntesection at each waypoint to extract lidar scans
        [pts, isOccupied] = rayIntersection(map3D, [gTruth(i,:) rot], dirAngles, 800);
        pts = cast(pts, 'like', tformPt.Location);
        ptsvalid = pts(isOccupied==1,:);
        orgPts = ptsvalid(ptsvalid(:,1)~=0 | ptsvalid(:,2)~=0 | ptsvalid(:,3)~=0,:);
        
        % Store original map for comparision
        mapPoints = [mapPoints;orgPts];
        
        % Transform lidar scan such that the waypoint will be the origin
        % for it
        tempPt = pointCloud(orgPts);
        tform = rigidtform3d([0 0 0],[-gTruth(i,1) -gTruth(i,2) -altitude]);
        pClouds(i) = pctransform(tempPt, tform);
    end
    
    orgMap = pointCloud(mapPoints);
end

function [loopSubmapIds,loopScores] = helperEstimateLoopCandidates(subMaps, currentScan, ...
    subMapPoses, recentPose, nsbmp, subMapThresh, loopClosureSearchRadius, loopClosureThreshold)
% The function matches current scan with submaps and tries to find the
% possible matches

    loopClosureMaxAttempts = 3;   % Maximum number of submaps checked
    maxLoopClosureCandidates = 3; % Maximum number of matches to send back

    %% Initialize variables
    loopSubmapIds = []; % Variable to hold the candidate submap IDs
    loopScores = [];    % Variable to hold the score of the matched submaps

    %% Find the submaps to be considered based on distance and Submap ID
    % Compute the centers of all Submaps
    centers = subMapPoses(1:nsbmp-subMapThresh,1:3);
    distanceToCenters = vecnorm(centers - recentPose(1:3),2,2);
    idx=find(distanceToCenters < loopClosureSearchRadius);
    centerCandidates = [distanceToCenters(idx) idx];

    % If there are submaps to be considered, sort them by distance from the
    % current scan
    if ~isempty(centerCandidates)
        % Sort them based on the distance from the current scan
        centerCandidates = sortrows(centerCandidates);

        % Return only the minimum number of loop candidates to be returned
        N = min(loopClosureMaxAttempts, size(centerCandidates,1));
        nearbySubmapIDs = centerCandidates(1:N, 2)';
    else
        nearbySubmapIDs = [];
    end
    %% Match the current scan with the candidate submaps
    newLoopCandidates = zeros(0,1); % Loop candidates
    newLoopCandidateScores = zeros(0,1); % Loop candidate scores
    count = 0; % Number of loop closure matches found

    % If there are submaps to consider
    if ~isempty(nearbySubmapIDs)
        % For every candidate submap
        for k = 1:length(nearbySubmapIDs)
            submapId = nearbySubmapIDs(k);
            [~,~,rmse] = pcregistericp(currentScan, subMaps{submapId}, 'MaxIterations',50);

            % Accept submap only if it meets the score threshold
            if rmse < loopClosureThreshold
                count = count + 1;
                % Keep track of matched Submaps and their scores
                newLoopCandidates(count) = submapId;
                newLoopCandidateScores(count) = rmse;
            end

        end

        % If there are candidates to consider
        if ~isempty(newLoopCandidates)
            % Sort them by their scores in descending order
            [~,ord] = sort(newLoopCandidateScores);
            % Return the required number of submaps matched, and their scores
            loopSubmapIds = newLoopCandidates(ord(1:min(count,maxLoopClosureCandidates)));
            loopScores = newLoopCandidateScores(ord(1:min(count,maxLoopClosureCandidates)));
        end
    end
end
