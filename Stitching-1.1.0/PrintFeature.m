function PrintFeature( input_A, input_B, all, back, outPath)
%PrintFeature print the feature correspondence.
%   Detailed explanation goes here
    fileListA = dir(input_A);
    fileListA = fileListA(3:length(fileListA));
    fileListB = dir(input_B);
    fileListB = fileListB(3:length(fileListB));
    nFrames = length(fileListA);
    if nFrames ~= size(all, 1)
        error('frame number not match');
    end
    for frameIndex = 1:nFrames
        IA = imread([input_A fileListA(frameIndex).name]);
        IB = imread([input_B fileListB(frameIndex).name]);
        pa_all = squeeze(all(frameIndex, :, 1:2));
        pb_all = squeeze(all(frameIndex, :, 3:4));
        pa_back = squeeze(back(frameIndex, :, 1:2));
        pb_back = squeeze(back(frameIndex, :, 3:4));

        IA = DrawFeature(IA, pa_all, 'red', 'o');
        IB = DrawFeature(IB, pb_all, 'red', 'o');

        IA = DrawFeature(IA, pa_back, 'green', 's');
        IB = DrawFeature(IB, pb_back, 'green', 's');
        if ~exist(outPath, 'dir')
            mkdir(outPath);
        end
        imwrite(IA, [outPath '/A' int2str(frameIndex) '.jpg']);
        imwrite(IB, [outPath '/B' int2str(frameIndex) '.jpg']);
    end
end

