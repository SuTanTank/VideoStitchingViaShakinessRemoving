function [ homos ] = NewWarping( pa_, pb_, H, W, qH, qW, lambda)
% A warpped up wersion of as-similar-as-possible warping, with pre-warp. 
% set PREWARP to false if the output is badly distorted. 
    PREWARP = true;
    nP = length(pa_);
    if length(pb_) ~= nP
        error('Points Numbers met Matching!');
    end
    if PREWARP
        [preH, ~] = ransacfithomography(pa_', pb_', 0.001);
        whilecount = 0;
        while isnan(preH(1, 1)) && whilecount < 100
            [preH, ~] = ransacfithomography(pa_', pb_', 0.001);
            whilecount = whilecount + 1;
        end
        if whilecount == 100
            error('?') ;
        end
    else
        preH = eye(3);
        pa = pa_;  
    end
    pbWarp = preH \ [pb_' ; ones(1, nP)];
    pbWarp(1, :) = pbWarp(1, :) ./ pbWarp(3, :);
    pbWarp(2, :) = pbWarp(2, :) ./ pbWarp(3, :);
    pbWarp = pbWarp(1:2, :)';    
    
    diff = sum((pa_ - pbWarp) .* (pa_ - pbWarp), 2) < 1000;    
    valid = pbWarp(:, 1) > 0 & pbWarp(:, 1) < W & pbWarp(:, 2) > 0 & pbWarp(:, 2) < H;
    pa = pa_(valid & diff, :);
    pbWarp = pbWarp(valid & diff, :);
    
    asap = AsSimilarAsPossibleWarping(H, W, qW, qH, lambda);
    asap.SetControlPts(pa, pbWarp);
    asap.Solve();
    
% -----DEBUG-SCRIPT------
% use it to check the warping result
%     e = asap.CalcError();
    
%     grid = asap.Warp(ones(H, W, 3) * 255, 400);
%     imshow(grid);
    
%     if e > 10
%         disp('?') ;
%     end

    homos2 = asap.CalcHomos();
    homos = homos2;
    for row = 1:H/qH
        for col = 1:W/qW
            tempH = preH * squeeze(homos2(row, col, :, :));
            homos(row, col, :, :) = tempH ./ tempH(3, 3);            
        end
    end
end

