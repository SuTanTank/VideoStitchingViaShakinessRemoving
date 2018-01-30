function L = coding_length(W, epsilon2, totalSampleCount, affine);

% coding_length.m
%
%   Estimates the number of bits needed to code a set of vectors upto
%   distortion epsilon2, based on a coding length derived from the optimal
%   Gaussian rate-distortion.
%
%   For details, see Section 3 and Appendices I and II of Ma et. al.
%   Segmentation of Multivariate Mixed Data via Lossy Coding and
%   Compression. PAMI 2007.
%
% Inputs:
%   W                   - [ w1 | w2 | ... | wm ], the matrix of data vectors (columns).
%   epsilon2            - the allowable squared distortion in representing the w_i
%                          E[ \| w_i - hat(w_i) \|^2 ] = epsilon2
%   totalSampleCount    - total number of samples in all groups (not just
%                          this one).
%   affine              - boolean, if true we use the coding length derived
%                          for nonzero-mean data.
%
% Outputs:
%   L - estimate of the number of bits needed to encode the
%        vectors.
%
% Dependencies:
%   none
%
% Feb. '06  John Wright -- jnwright@uiuc.edu
% Sep. '07  Shankar Rao -- srrao@uiuc.edu

% Copyright 2007, University of Illinois. All rights reserved.

[n, m] = size(W);
if affine
    mu_hat = mean(W,2);
    W = W - repmat(mu_hat, 1, m);
    L = (n/2)*log2(1+sum(mu_hat.^2)/epsilon2) ...
        + ((m+n)/2)*sum(log2(1+(n/(epsilon2*m))*svd(W, 0).^2)) ...
        - m*log2(m/totalSampleCount);
else
    L = ((m+n)/2)*sum(log2(1+(n/(epsilon2*m))*svd(W, 0).^2)) - m*log2(m/totalSampleCount);
end
