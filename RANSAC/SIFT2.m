function [source_features,target_features, score, source_features_o, target_features_o] = SIFT2(I1,I2, validMask)

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
    source_features_o = kp1(1:2,matches(1,:))';
    target_features_o = kp2(1:2,matches(2,:))';
    
    % give the voting index row first 20 * 20
    
    [H, W, ~] = size(I1);
    qH = H / 10;
    qW = W / 10;
    
    index = floor(source_features_o(:, 1) / qW) * 10 + (floor(source_features_o(:, 2) / qH)) + 1;
    
    [~,inliners, score] = EstimateHomographyByRANSAC2(source_features_o',target_features_o', 0.001, index, validMask);
%     [~,inliners] = EstimateHomographyByRANSAC(source_features_o',target_features_o', 0.001);
    source_features = source_features_o(inliners,:);
    target_features = target_features_o(inliners,:);

    
    
    
end