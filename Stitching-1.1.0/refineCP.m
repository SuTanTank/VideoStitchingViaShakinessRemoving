function [ CP, ppf ] = refineCP( CP_, ppf_, backListA, backListB, trackA, trackB)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    nFrame = trackA.nFrame;
    CP = zeros(size(CP_));
    ppf = zeros(nFrame, 1);
    videoH = trackA.videoHeight;
    videoW = trackA.videoWidth;
    for frameIndex = 1:nFrame
        fprintf('%5d', frameIndex);
        if mod(frameIndex, 20) == 0
            fprintf('\n') ;
        end
        maskA = trackA.getBackMask(frameIndex, backListA);
        maskB = trackB.getBackMask(frameIndex, backListB);
%         figure(1);imshow(maskA);
%         figure(2);imshow(maskB);
%         if frameIndex == 450
%             imwrite(maskA, 'D:/450-maska.jpg') ;
%             imwrite(maskB, 'D:/450-maskb.jpg') ;
%         end
        pa = squeeze(CP_(frameIndex, 1:ppf_(frameIndex), 1:2));
        pb = squeeze(CP_(frameIndex, 1:ppf_(frameIndex), 3:4));
        valid = pa(:, 1)>0.5 & pa(:, 1)<videoW-0.5 & pb(:, 1)> 0.5 & pb(:, 1)<=videoW-0.5 ...
            & pa(:, 2)>0.5 & pa(:, 2)<=videoH-0.5 & pb(:, 2)> 0.5 & pb(:, 2)<=videoH-0.5;
        pa = pa(valid, :);
        pb = pb(valid, :);
        cpIndex = maskA(round(pa(:, 1) - 1) * videoH + round(pa(:, 2))) > 0.5 ...
            & maskB(round(pb(:, 1) - 1) * videoH + round(pb(:, 2))) > 0.5;
        ppf(frameIndex) = sum(cpIndex);
        CP(frameIndex, 1:ppf(frameIndex), 1:2) = pa(cpIndex, :);
        CP(frameIndex, 1:ppf(frameIndex), 3:4) = pb(cpIndex, :);
    end
end

