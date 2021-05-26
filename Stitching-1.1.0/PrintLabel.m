function PrintLabel( input, tracks, savePath)
    fileList = dir(input);
    fileList = fileList(3:length(fileList));
    nFrames = length(fileList);
    if ~exist(savePath, 'dir')
        mkdir(savePath);
    end
    if nFrames ~= tracks.nFrame
        error('frame number not match');
    end
    colors = {'black', 'red', 'green', 'magenta', 'blue', 'cyan', 'yellow', 'blue', 'cyan', 'yellow'};
    shapes = {'x', 'o', '+', '*', 's', 'o', '+', '*', 's', 'o'};
    label_num = 9;
    for frameIndex = 1:nFrames
        fileName = fileList(frameIndex).name;
        frame = imread([input fileName]);
        hasPoint = tracks.points(:, frameIndex, 1) ~= 0 | tracks.points(:, frameIndex, 2) ~= 0;
        points = squeeze(tracks.points(hasPoint, frameIndex, :));
        labels = squeeze(tracks.labels(hasPoint, frameIndex, :));
        I1 = frame;
        I2 = frame;
        for labelIndex = 1:label_num
            pos = squeeze(points(labels(:, 1) == labelIndex, :));
%             I1 = insertMarker(I1, pos, shapes{labelIndex+1}, 'color', colors{labelIndex+1});
            I1 = DrawFeature(I1, pos, colors{labelIndex+1}, 'o');
            pos = squeeze(points(labels(:, 2) == labelIndex, :));
%             I2 = insertMarker(I2, pos, shapes{labelIndex+1}, 'color', colors{labelIndex+1});
            I2 = DrawFeature(I2, pos, colors{labelIndex+1}, 'o');
        end
%         figure(1); imshow(I1);
%         figure(2); imshow(I2);
        imwrite(I1, [savePath int2str(frameIndex) '_1.jpg']);
        imwrite(I2, [savePath int2str(frameIndex) '_2.jpg']);
    end

end

