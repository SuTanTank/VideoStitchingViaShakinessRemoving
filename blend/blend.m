function [im] = blend(im1, im2, c1, c2)

[height, width, ~] = size(im1);

alpha1 = zeros(size(im1));
alpha2 = zeros(size(im2));

im1 = im2double(im1);
im2 = im2double(im2);


im = (alpha2 .* im2 + alpha1 .* im1) ./ (alpah1 + alpha2);