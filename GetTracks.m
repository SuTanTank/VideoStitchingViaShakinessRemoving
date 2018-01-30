function tracks = GetTracks( input, meshSize, demand)
%GetTracks Compute tracks by KLT
%   Use KLT to track evenly fistributed track points
%   input: the path to images
%   meshSize: the meshSize of video stitching
    fileList = dir(input);
    fileList = fileList(3:length(fileList));
    nFrames = length(fileList); 
    
    tracks = TrackLib(nFrames);
    
    tracker = vision.PointTracker('MaxBidirectionalError', 1);
    fileName = fileList(1).name;
    frame = imread([input fileName]);
    [H, W, ~] = size(frame);
    tracks.videoWidth = W;
    tracks.videoHeight = H;
    livePoints = getMorePoints(frame, meshSize, 0, [], demand);
%     livePoints = filtermask(frame, livePoints);
    initialize(tracker, livePoints, frame);
    tracks.addPoints(livePoints, 1);
    for frameIndex = 2:length(fileList)
        fprintf('%5d', frameIndex);
        if mod(frameIndex, 20) == 0
            fprintf('\n') ;
        end        
        fileName = fileList(frameIndex).name;
        frame = imread([input fileName]);
        
        [livePoints, validity] = step(tracker, frame);
        age = true(size(validity));
        age(tracks.len(tracks.live) == 100) = false;
%         fprintf('=> %d\t%d\n', size(tracks.live, 2), size(validity, 1));
        if size(tracks.live, 2) ~= size(validity, 1)
            disp('?') ;
        end
        tracks.endPoints(validity & age, frameIndex);
        tracks.updatePoints(livePoints(validity & age, :), frameIndex);
        
        % end too old tracks 
        
        
        morePoints = getMorePoints(frame, meshSize, length(tracks.live), livePoints(validity == true, :), demand);
%         morePoints = filtermask(frame, morePoints);
        tracks.addPoints(morePoints, frameIndex);
        livePoints = [livePoints(validity & age, :); morePoints];
        setPoints(tracker, livePoints);
        marked = insertMarker(frame, livePoints);
        imshow(marked);
        
    end
    tracks.endPoints(false(length(tracks.live), 1), length(fileList) + 1);
end

function pointsMore = getMorePoints(frame, meshSize, nP, oldpoints, demand)
    demand = demand / (meshSize * meshSize);
    votes = zeros(meshSize);
    [H, W, ~] = size(frame);
    threshold = 0.5;
    if nP > 0
        votes = getVotes(frame, meshSize, oldpoints);
    end
    points = [];
    
    for row = 1:meshSize
        for col = 1:meshSize
            if votes(row, col) < demand * 0.7
                nMore = demand - votes(row, col);
                roi = [1 + (col - 1) * W / meshSize, 1 + (row - 1) * H / meshSize, W / meshSize - 1, H / meshSize - 1];                
                pNew = detectMinEigenFeatures(rgb2gray(frame), 'ROI', roi, 'MinQuality', threshold); 
                while (size(pNew, 1) < nMore) && threshold > 0.1
                    threshold = threshold - 0.1; 
                    threshold = max(threshold, 0);
                    pNew = detectMinEigenFeatures(rgb2gray(frame), 'ROI', roi, 'MinQuality', threshold); 
                end
                if nMore < size(pNew, 1)
                    ordering = randperm(length(pNew));
                    pNew = pNew(ordering);
                    pNew = pNew(1:nMore);
                end
                points = [points; pNew.Location];
            end
        end
    end
    pointsMore = points;
    
end

function votes = getVotes(frame, meshSize, points)
    [H, W, ~] = size(frame);    
    qH = H / meshSize;
    qW = W / meshSize;
    index = floor(points(:, 1) / qW) * meshSize + (floor(points(:, 2) / qH)) + 1; 
    voting = histcounts([index; 1; meshSize*meshSize], meshSize*meshSize);  
    voting(1) = voting(1) - 1;
    voting(meshSize * meshSize) = voting(meshSize*meshSize) - 1;    
    votes = reshape(voting, [meshSize meshSize]);
%     votes = votes';    
end

function points = filtermask(frame, points_)
    mask = frame(:, :, 1) < 20 & frame(:, :, 2) < 20 & frame(:, :, 3) < 20;
%     mask_ = ~mask;
    mask = imgaussfilt(double(mask), 50);
%     mask_ = imgaussfilt(double(mask_), 10);
%     mask = double(mask_ > mask);
    videoH = size(frame, 1);
%     imshow(mask * 100);
    mask(mask > 0.2) = 1;
    mask(mask ~= 1) = 0;
%     imshow(mask);
    valid = mask(round(points_(:, 1) - 1) * videoH + round(points_(:, 2))) == 0;
    points = points_(valid, :);
end
