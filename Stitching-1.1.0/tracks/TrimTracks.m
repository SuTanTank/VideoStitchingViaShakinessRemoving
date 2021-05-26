function [trackA, trackB] = TrimTracks(trackA, trackB)
nFrame = min(trackA.nFrame, trackB.nFrame);

goodA = trackA.head <= nFrame; 
trackA.points = trackA.points(goodA, 1:nFrame, :);
trackA.head = trackA.head(goodA);
trackA.tail = trackA.tail(goodA);
trackA.tail = min(trackA.tail, nFrame);
trackA.len = trackA.tail - trackA.head + 1;
trackA.nFrame = nFrame;
trackA.nTrack = length(trackA.head);

goodB = trackB.head <= nFrame; 
trackB.points = trackB.points(goodB, 1:nFrame, :);
trackB.head = trackB.head(goodB);
trackB.tail = trackB.tail(goodB);
trackB.tail = min(trackB.tail, nFrame);
trackB.len = trackB.tail - trackB.head + 1;
trackB.nFrame = nFrame;
trackB.nTrack = length(trackB.head);

end

