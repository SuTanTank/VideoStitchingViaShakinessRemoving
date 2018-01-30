function [source_features,target_features] = SIFT(I1,I2)
    [ kp1,ds1 ] = vl_sift(single(rgb2gray(I1)),'PeakThresh', 1,'edgethresh',500);
    [ kp2,ds2 ] = vl_sift(single(rgb2gray(I2)),'PeakThresh', 1,'edgethresh',500);
    
    i = 1;
    while i < length(kp1)
        if kp1(1, i) - kp1(1, i+1) < 0.1 && kp1(2, i) - kp1(2, i+1) < 0.1
            kp1(:, i) = [];
            ds1(:, i) = [];
        else
            i = i + 1;
        end            
    end
    i = 1;
    while i < length(kp2)
        if kp2(1, i) - kp2(1, i+1) < 0.1 && kp2(2, i) - kp2(2, i+1) < 0.1
            kp2(:, i) = [];
            ds2(:, i) = [];
        else
            i = i + 1;
        end            
    end
    
    matches   = vl_ubcmatch(ds1,ds2);
    ordering = randperm(length(matches));
    matches = matches(:, ordering);
    source_features = kp1(1:2,matches(1,:))';
    target_features = kp2(1:2,matches(2,:))';
    
%     [~,inliners] = EstimateHomographyByRANSAC(source_features',target_features', 0.001);
%     source_features = source_features(inliners,:);
%     target_features = target_features(inliners,:);    
end