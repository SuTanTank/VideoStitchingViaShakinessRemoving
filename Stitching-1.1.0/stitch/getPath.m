function [ path ] = getPath( input, MeshSize, tracks, backList)
%GETPATH Compute the bundled camera path
%   track is the set of trajectories of the input video, backList indicate
%   the background trajectories. 
    disp(['Computing camera path of ' input]);
    if input(length(input)) ~= '/'
        input = [input '/'] ;
    end
    
    window_size = 40;
    tracks.setWindowSize(window_size);    
    
    fileList = dir(input);
    fileList = fileList(3:length(fileList));
    nFrames = tracks.nFrame;
    if nFrames < 2
        error('Wrong inputs') ;
    end
    path = zeros(nFrames, MeshSize, MeshSize, 3, 3);
    for row = 1:MeshSize
        for col = 1:MeshSize
            path(1, row, col, :, :) = eye(3);
        end
    end
    fprintf('%5d', 1);
    for frameIndex = 2:nFrames
        fprintf('%5d', frameIndex);
        if mod(frameIndex, 20) == 0
            fprintf('\n') ;
        end
        fileName = fileList(frameIndex - 1).name;
        I1 = imread([input fileName]);
        fileName = fileList(frameIndex).name;
        I2 = imread([input fileName]);
        [H, W, ~] = size(I1);
        quadH = H / MeshSize;
        quadW = W / MeshSize;
        asaplambda = 1;
        [I1_features,I2_features] = tracks.getGoodF(frameIndex - 1, backList);
        homos = NewWarping(I1_features, I2_features, H, W, quadH, quadW, asaplambda);
        if length(I1_features) < 4
            error('not enough matched features');
        end
        for i = 1:MeshSize
            for j = 1:MeshSize
                path(frameIndex, i, j, :, :) = squeeze(homos(i, j, :, :)) * squeeze(path(frameIndex - 1, i, j, :, :));
                path(frameIndex, i, j, :, :) = squeeze(path(frameIndex, i, j, :, :)) / path(frameIndex, i, j, 3, 3);
            end
        end            
    end
end

