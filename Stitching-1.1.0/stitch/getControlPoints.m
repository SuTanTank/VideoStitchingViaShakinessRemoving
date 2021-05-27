function [CP, ppf] = getControlPoints( input_A, input_B, maxppf, ransac)
    disp('Detecting and Matching SIFT features...');
    fileListA = dir(input_A);
    fileListA = fileListA(3:length(fileListA));
    fileListB = dir(input_B);
    fileListB = fileListB(3:length(fileListB));
    nFrames = min(length(fileListA), length(fileListB));
    CP = zeros(nFrames, maxppf, 4);
    ppf = zeros(nFrames, 1);
    trackerA = vision.PointTracker('MaxBidirectionalError', 1);
    trackerB = vision.PointTracker('MaxBidirectionalError', 1);
    for frameIndex = 1:nFrames
        fprintf('%5d', frameIndex);
        if mod(frameIndex, 20) == 0
            fprintf('\n') ;
        end

        fileNameA = fileListA(frameIndex).name;
        fileNameB = fileListB(frameIndex).name;
        IA = imread([input_A fileNameA]);
        IB = imread([input_B fileNameB]);

        [H, W, ~] = size(IA);

        if frameIndex > 1
            setPoints(trackerA, trackA);
            [trackAcont, validityA] = step(trackerA, IA);
            setPoints(trackerB, trackB);
            [trackBcont, validityB] = step(trackerB, IB);
            trackAcont = trackAcont(validityA & validityB, :);
            trackBcont = trackBcont(validityA & validityB, :);
        end
        if ransac
            [trackAsurf, trackBsurf] = SURF(IA, IB);
        else
            [trackAsurf, trackBsurf] = SURF2(IA, IB);
        end

        trackA = trackAsurf;
        trackB = trackBsurf;
        
        valid = trackA(:, 1) > 0 & trackA(:, 1) < W & trackA(:, 2) > 0 & trackA(:, 2) < H ...
            & trackB(:, 1) > 0 & trackB(:, 1) < W & trackB(:, 2) > 0 & trackB(:, 2) < H;
        
        trackA = trackA(valid, :);
        trackB = trackB(valid, :);

        valid = filtermask(IA, trackA);
        trackA = trackA(valid, :);
        trackB = trackB(valid, :);
        valid = filtermask(IB, trackB);
        trackA = trackA(valid, :);
        trackB = trackB(valid, :);

        if frameIndex == 1
            initialize(trackerA, trackA, IA);
            initialize(trackerB, trackB, IB);
        else
            if length(trackA) + length(trackAcont) > maxppf
                ordering = randperm(length(trackAcont));
                trackAcont = trackAcont(ordering(1:maxppf - length(trackA)), :);
                trackBcont = trackBcont(ordering(1:maxppf - length(trackA)), :);
            end
            trackA = [trackA ; trackAcont];
            trackB = [trackB ; trackBcont];
        end
        
        IA = insertMarker(IA, trackA, 'o', 'color', 'red');
        IB = insertMarker(IB, trackB, 's', 'color', 'yellow');
        figure(1);
        imshow(IA);
        figure(2);
        imshow(IB);

        if length(trackA) > maxppf
            ppf(frameIndex) = maxppf;
            CP(frameIndex, :, 1:2) = trackA(1:maxppf, :);
            CP(frameIndex, :, 3:4) = trackB(1:maxppf, :);
        else
            ppf(frameIndex) = length(trackA);
            CP(frameIndex, 1:ppf(frameIndex), 1:2) = trackA;
            CP(frameIndex, 1:ppf(frameIndex), 3:4) = trackB;
        end
%         CP(frameIndex, :, :) = [featuresA featuresB; zeros(maxppf - ppf(frameIndex), 4)];

    end
end

function valid = filtermask(frame, points_)
    [H, W, ~] = size(frame);
    valid = points_(:, 1) > 0 & points_(:, 1) < W & points_(:, 2) > 0 & points_(:, 2) < H;
    points_ = points_(valid, :);
    mask = frame(:, :, 1) < 20 & frame(:, :, 2) < 20 & frame(:, :, 3) < 20;
    mask = imgaussfilt(double(mask), 50);
    videoH = size(frame, 1);
    mask(mask > 0.2) = 1;
    mask(mask ~= 1) = 0;
    valid = mask(round(points_(:, 1) - 1) * videoH + round(points_(:, 2))) == 0;
end


