function codingLength = total_coding_length(data, labels, epsilon, affine)

% total_coding_length
%
%   compute the total number of bits needed to code a set of data vectors
%   given a grouping, and an error tolerance.
%
% Inputs:
%   data - a matrix where each column is a data vector.
%   labels - a list of integer labels that partitions the data into groups
%   epsilon - the error tolerance
%   affine - if true, compute using the affine coding length function.
%               otherwise use the linear coding length function
%
% Outputs:
%   codingLength - the total number of bits needed to code this set of data
%
% Dependencies:
%   coding_length
%
% Sep. '07  Shankar Rao -- srrao@uiuc.edu

% Copyright 2007, University of Illinois. All rights reserved.

codingLength = 0;
groupCount = max(labels);
sampleCount = length(labels);
% Make sure outlier indices are counted too
% modified 10-24-07
setOfIndices = unique(labels);
for groupIndex = setOfIndices(:)',
    codingLength = codingLength + coding_length(data(:, labels == groupIndex), epsilon^2, sampleCount, affine);
end