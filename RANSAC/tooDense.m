function [ res ] = tooDense( x, threshold)
%tooDense Summary of this function goes here
%   Examine whether 4 samples are too close to each other
    x1 = x(1:2, :);
    x2 = x(4:5, :);
    score = min(max(x1(1, :)) - min(x1(1, :)) + max(x1(2, :)) - min(x1(2, :)), ...
        max(x2(1, :)) - min(x2(1, :)) + max(x2(2, :)) - min(x2(2, :)));
    if score < threshold
        res = 1;
    else
        res = 0;
    end

end

