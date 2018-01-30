function backList = refineTrack( track, good )
%GETREST Refinement of trajectories
    good = [ones(track.nWindow, 1) * (-1) good];    
    backList = zeros(track.nTrack, 1);
    trackLength = zeros(track.nTrack, track.nWindow);
    trackWindowGood = zeros(track.nTrack, track.nWindow);
    for frameIndex = 1:track.nFrame
        windowIndex = ceil(frameIndex / (track.wSize / 2));        
        if windowIndex < track.nWindow + 1
            hasLabel = track.labels(:, frameIndex, 2) > 0;
            hasGoodLabel = good(windowIndex, track.labels(:, frameIndex, 2) + 1) == 1;
            trackLength(hasLabel, windowIndex) = 1;
            trackWindowGood(hasGoodLabel, windowIndex) = 1;
        end
        if windowIndex > 1
            hasLabel = track.labels(:, frameIndex, 1) > 0;
            hasGoodLabel = good(windowIndex - 1, track.labels(:, frameIndex, 1) + 1) == 1;
            trackLength(hasLabel, windowIndex - 1) = 1;
            trackWindowGood(hasGoodLabel, windowIndex - 1) = 1;
        end
    end
    backList((sum(trackWindowGood, 2) - ceil( 0.7 * sum(trackLength, 2)) >= 0) & (sum(trackLength, 2) > 0)) = 1;
    backList((sum(trackWindowGood, 2) - ceil( 0.7 * sum(trackLength, 2)) < 0) & (sum(trackLength, 2) > 0)) = -1;
    backList(sum(trackWindowGood, 2) - sum(trackLength, 2) > 0) = -10;
    
end

