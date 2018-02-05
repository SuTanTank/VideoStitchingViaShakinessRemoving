% Video Stitching for videos captured by handheld device
% Written by Tan SU
% contact: sutank922@gmail.com
addpath('mesh');
addpath('RANSAC');
addpath('stitch');
addpath('blend');
addpath('tracks');
addpath('tracks/helpers');
addpath('graph');
addpath('peter');
%% Set up vlfeat
poolobj = gcp('nocreate');
delete(poolobj);
cd vlfeat-0.9.19/toolbox;
feval('vl_setup');
cd ../..;

%% params
data = '../case-cuhk_lib/';
input_A = 'left/';
input_B = 'right/';
% ---------------
PointsPerFrame = 500; % # of feature correspondences in a frame
TracksPerFrame = 600; % # of trajectories in a frame
TrackWindowSize = 40; % the window size for motion segmentation
% ---------------
MeshSize = 8; % The mesh size of bundled camera path, 5 - 10 is OK
MaxIte = 15; % Number of iterations of the optimization scheme, 10 - 15 is OK
Smoothness = 1; % adjust how stable the output is, 0.5 - 3 is OK
Cropping = 1; % adjust how similar the result to the original video, usually set to 1;
Stitchness = 20; % adjust the weight of stitching term, 10 - 30 is OK
% ---------------
OutputPadding = 500; % the padding around the video
OutputPath = 'res_demo'; % the directory to store the output frames, auto create it if not exist

%% intermediate output control
PRINT_LABEL = true;
PRINT_BACKGROUND = true;
PRINT_FEATURES = true;
PRINT_GRID = true;


%% Track by KLT
tic;
if ~exist([data 'tracks' int2str(TracksPerFrame) '.mat'], 'file')
    trackA = GetTracks([data input_A], 10, TracksPerFrame); % 10 = tracks evenly distributed in 10*10 grids
    trackB = GetTracks([data input_B], 10, TracksPerFrame);
    save([data 'tracks' int2str(TracksPerFrame) '.mat'], 'trackA', 'trackB');

else
    load([data 'tracks' int2str(TracksPerFrame) '.mat']);
end
toc;

%% Clustering trajectories
tic;
if ~exist([data 'tracks_' int2str(TrackWindowSize) '_' int2str(TracksPerFrame) '.mat'], 'file')
    trackA.addLabel(40);
    trackB.addLabel(40);
    if trackA.nLabel > trackB.nLabel
        trackB.nLabel = trackA.nLabel;
    else
        trackA.nLabel = trackB.nLabel;
    end

    save([data 'tracks_' int2str(TrackWindowSize) '_' int2str(TracksPerFrame) '.mat'], 'trackA', 'trackB');
else
    load([data 'tracks_' int2str(TrackWindowSize) '_' int2str(TracksPerFrame) '.mat']);
end
toc;

%% print the label
if PRINT_LABEL
    PrintLabel([data input_A], trackA, [data '/left_label2/']);
    PrintLabel([data input_B], trackB, [data '/right_label2/']);
end

%% Matching SIFT in every frame pair
tic;
if ~exist([data 'ControlPoints' int2str(PointsPerFrame) '.mat'], 'file')
    [CP, ppf] = getControlPoints([data input_A], [data input_B], 500);
    save([data 'ControlPoints' int2str(PointsPerFrame) '.mat'], 'CP', 'ppf');
else
    load([data 'ControlPoints' int2str(PointsPerFrame) '.mat']);
end
toc;
%% Get common background
tic;
if ~exist([data 'graph' int2str(TracksPerFrame) '.mat'], 'file')
    alpha = 0.01;
    beta = 0.01; % these two parameters are used in setting the graph edges
    maxNlabel = max(trackA.nLabel, trackB.nLabel);
    trackA.nLabel = maxNlabel; trackB.nLabel = maxNlabel;
    [path, graph, goodA, goodB] = GetGraph(trackA, trackB, CP, ppf, alpha, beta);
    backListA = refineTrack(trackA, goodA);
    backListB = refineTrack(trackB, goodB);
    save([data 'graph' int2str(TracksPerFrame) '.mat'], 'path', 'graph', 'backListA', 'backListB');    
else
    load([data 'graph' int2str(TracksPerFrame) '.mat']);
end
if PRINT_BACKGROUND
    PrintBackground([data '/left/'], trackA, backListA, [data '/left_back/']);
    PrintBackground([data '/right/'], trackB, backListB, [data '/right_back/']);
end
toc;
%% Compute original camera path (by As-similar-as-possible Warping)
% the rigidity can be controlled by setting asaplambda inside getPath.m
tic;
if ~exist([data 'Path' int2str(MeshSize) '.mat'], 'file')
    [pathA] = getPath([data input_A], MeshSize, trackA, backListA);
    [pathB] = getPath([data input_B], MeshSize, trackB, backListB);
    save([data 'Path' int2str(MeshSize) '.mat'], 'pathA', 'pathB');
else
    load([data 'Path' int2str(MeshSize) '.mat']);
end
toc;
tic;
if ~exist([data 'ControlPoints' int2str(PointsPerFrame) '_refine.mat'], 'file')
    [CP_refine, ppf_refine] = refineCP(CP, ppf, backListA, backListB, trackA, trackB);
    save([data 'ControlPoints' int2str(PointsPerFrame) '_refine.mat'], 'CP_refine', 'ppf_refine');
else
    load([data 'ControlPoints' int2str(PointsPerFrame) '_refine.mat']);
end
toc;
if PRINT_FEATURES
    PrintFeature([data input_A], [data input_B], CP, CP_refine, [data '/features/']);
end


%% Optimize the paths
tic;
stitcher = VideoStitch2([data input_A], [data input_B], pathA, pathB, CP, ppf, Smoothness, Cropping, Stitchness);
stitcher.init();
SECOND_ROUND = 5; 
% use it to perform a 2nd phase optimization for more stable output at non-overlapping region
% set to 5 - 10 is OK, to turn off the 2nd phase optimization, set it to 0.
% WARNING: should be larger than MaxIte
stitcher.optPath(MaxIte, SECOND_ROUND);

% save([data 'stitcher_no.mat'], 'stitcher');
stitcher.render([data OutputPath], OutputPadding);
if PRINT_GRID
    stitcher.renderGrid([data '/res_grid/'], OutputPadding);
end
toc;
%% Evaluation Score
padding = 0;
CP_ = stitcher.getStitchedCP(padding + 1, stitcher.nFrames - padding);
save([data 'stitchedCP.mat'], 'CP_');