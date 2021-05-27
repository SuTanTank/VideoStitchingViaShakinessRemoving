function [source_features,target_features] = SURF2(I1,I2)

    [H, W, ~] = size(I1);
    grayI1 = rgb2gray(I1);
    grayI2 = rgb2gray(I2);

%     [f1, vpts1] = getSURFFeatures(grayI1);
%     [f2, vpts2] = getSURFFeatures(grayI2);    
    threshold = 10;
    points1 = detectSURFFeatures(grayI1, 'MetricThreshold', threshold);
    points2 = detectSURFFeatures(grayI2, 'MetricThreshold', threshold);
    [f1, vpts1] = extractFeatures(grayI1, points1);  
    [f2, vpts2] = extractFeatures(grayI2, points2);  
    
    index_pairs = matchFeatures(f1, f2) ;
    matched_pts1 = vpts1(index_pairs(:, 1), :);
    matched_pts2 = vpts2(index_pairs(:, 2), :);

    [n,~] = size(matched_pts1);

    source_features = zeros(n,2);
    target_features = zeros(n,2);

    for i=1:n
        source_features(i,:) = matched_pts1(i, :).Location;
        target_features(i,:) = matched_pts2(i, :).Location;
    end   
    
    valid = source_features(:, 1) > 0 & source_features(:, 2) > 0 ...
        & source_features(:, 1) < W & source_features(:, 2) < H ...
        & target_features(:, 1) > 0 & target_features(:, 2) > 0 ...
        & target_features(:, 1) < W & target_features(:, 2) < H;
    source_features = source_features(valid, :);
    target_features = target_features(valid, :);
%     source_features = matched_pts1;
%     target_features = matched_pts2;
end

function [f, vpts] = getSURFFeatures(I)
    meshSize = 1;
    threshold = 100;
    [H, W] = size(I);
    f = [];
    vpts = [];
    for row = 1:meshSize
        for col = 1:meshSize            
            nMore = 1000;                
            roi = [1 + (col - 1) * W / meshSize, 1 + (row - 1) * H / meshSize, W / meshSize, H / meshSize];                
            pNew = detectSURFFeatures(I, 'ROI', roi, 'MetricThreshold', threshold); 
            while (size(pNew, 1) < nMore) && threshold > 10
                threshold = threshold - 10;                    
                pNew = detectSURFFeatures(I, 'ROI', roi, 'MetricThreshold', threshold); 
            end
            while nMore * 2 < size(pNew, 1) && threshold < 200
                threshold = threshold + 20;                    
                pNew = detectSURFFeatures(I, 'ROI', roi, 'MetricThreshold', threshold);                     
            end
            if nMore < size(pNew, 1)
                ordering = randperm(length(pNew));
                pNew = pNew(ordering);
                pNew = pNew(1:nMore);
            end               
            [fNew, vptsNew] = extractFeatures(I, pNew);  
            f = cat(1, f, fNew);
            vpts = cat(1, vpts, vptsNew.Location);
        end
    end
end