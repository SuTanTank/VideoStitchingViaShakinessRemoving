function [CP, ppf] = getControlPoints( input_A, input_B, maxppf )
    disp('Detecting and Matching SIFT features...');
    fileListA = dir(input_A);
    fileListA = fileListA(3:length(fileListA));
    fileListB = dir(input_B);
    fileListB = fileListB(3:length(fileListB));
    nFrames = length(fileListA);
    if nFrames ~= length(fileListB)
        error('wrong inputs');
    end
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
        
%         [~, ~, HH] = SURF(IA, IB);
        HH = eye(3);
%         IApre = imwarp(IA, projective2d(HH'), 'OutputView', imref2d(size(IA)));
        IApre = IA;
%         [trackApresift, trackBsift] = SIFT(IApre, IB);
        
        [trackApresurf, trackBsurf] = SURF2(IApre, IB);
%         trackerAB = vision.PointTracker('MaxBidirectionalError', 0.1);
%         trackApre = getMorePoints(IApre, 20, 2000);
%         initialize(trackerAB, trackApre, IApre);
%         [trackB, validity] = step(trackerAB, IB);
        
%         if length(trackB) > maxppf - length(trackBsurf)
%             ordering = randperm(length(trackB));
%             trackApre = trackApre(ordering(1:maxppf - length(trackBsurf)), :); 
%             trackB = trackB(ordering(1:maxppf - length(trackBsurf)), :); 
%             validity = validity(ordering(1:maxppf - length(trackBsurf)));
%         end
        
%         trackB = [trackB(validity, :); trackBsurf];
%         trackApre = [trackApre(validity, :); trackApresurf];

%         trackApre = [trackApresurf;trackApresift];
%         trackB = [trackBsurf;trackBsift];
        trackApre = trackApresurf;
        trackB = trackBsurf;
        trackA = HH \ [trackApre' ; ones(1, length(trackApre))];
        trackA(1, :) = trackA(1, :) ./ trackA(3, :);
        trackA(2, :) = trackA(2, :) ./ trackA(3, :);
        trackA = trackA(1:2, :);
        trackA = trackA';
        
        valid = trackA(:, 1) > 0 & trackA(:, 1) < W & trackA(:, 2) > 0 & trackA(:, 2) < H ...
            & trackB(:, 1) > 0 & trackB(:, 1) < W & trackB(:, 2) > 0 & trackB(:, 2) < H;
        trackA = trackA(valid, :);
        trackB = trackB(valid, :);
        
%         valid = filtermask(IA, trackA);
%         trackA = trackA(valid, :);
%         trackB = trackB(valid, :);
%         valid = filtermask(IB, trackB);
%         trackA = trackA(valid, :);
%         trackB = trackB(valid, :);
        
        if frameIndex == 1
            initialize(trackerA, trackA, IA); 
            initialize(trackerB, trackB, IB);
        else
            trackA = [trackA ; trackAcont];
            trackB = [trackB ; trackBcont];
        end
        if length(trackA) > maxppf
            ordering = randperm(length(trackA));
            trackA = trackA(ordering(1:maxppf), :); 
            trackB = trackB(ordering(1:maxppf), :); 
        end
        
%         IA = insertMarker(IA, trackA, 'o', 'color', 'red');
%         IB = insertMarker(IB, trackB, 's', 'color', 'yellow');
%         figure(1);
%         imshow(IA);
%         figure(2);
%         imshow(IB);
%         figure(3);
%         imshow(IApre);
        
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


