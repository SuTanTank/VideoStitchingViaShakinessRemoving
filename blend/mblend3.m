function [im] = mblend3(im1, im2, alpha)

[height, width, ~] = size(im1);

mask1 = im1;
mask2 = im2;

mask1 = mask1 < 1;
mask2 = mask2 < 1;

% r = 10;
% mask1 = dilate(mask1, r);
% mask2 = dilate(mask2, r);

x1 = mean(im1,3);
x1 = mean(x1, 1);
x1 = find(x1);
x1min = min(x1);
x1max = max(x1);

x2 = mean(im2,3);
x2 = mean(x2, 1);
x2 = find(x2);
x2min = min(x2);
x2max = max(x2);

if x1min < x2min
    left_x = x2min;
    right_x = x1max;
    alpha = sigmf(left_x:1:right_x,[0.02 (left_x + right_x)/2]);
    alpha = [zeros(1, left_x-1) alpha ones(1, width - right_x)];
else
    left_x = x1min;
    right_x = x2max;
    alpha = 1 - sigmf(left_x:1:right_x,[0.02 (left_x + right_x)/2]);
    alpha = [ones(1, left_x-1) alpha zeros(1, width - right_x)];
end
alpha = repmat(alpha, height, 1, 3);


im1(mask1) = im2(mask1);
im2(mask2) = im1(mask2);

im1 = im2double(im1);
im2 = im2double(im2);


im = alpha .* im2 + (1 - alpha) .* im1;

