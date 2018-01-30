function [distinctLabels labelIndices] = distinct_labels(labelMatrix)

% distinctLabels.m
%
%   Find the set of unique segmentations in a list of segmentations
%
% Inputs:
%   labelMatrix         - a matrix where each column corresponds to a
%                           labelling of data samples
%
% Outputs:
%   distinctLabels      - a matrix where each column corresponds to a 
%                           labelling of data samples, and no labelling is
%                           repeated.
%   labelIndices        - a vector where the ith element contains the index
%                         of the ith labeling in distinctLabels 
% Dependencies:
%   compare_labels
%
% Mar. '08  Shankar Rao -- srrao@uiuc.edu
%
% Copyright 2008, University of Illinois. All rights reserved.


if isempty(labelMatrix)
    distinctLabels = [];
    return;
end

trialCount = size(labelMatrix, 2);
labelIndices = zeros(1, trialCount);
distinctLabels = labelMatrix(:,1);
labelCount = 1;
labelIndices(1) = 1;

for trialIndex = 2:trialCount
    currentLabels = labelMatrix(:, trialIndex);
    foundMatch = false;
    for labelIndex = 1:labelCount
        err = compare_labels(currentLabels', distinctLabels(:, labelIndex)');
        if (err == 0)
            foundMatch = true;
            labelIndices(trialIndex) = labelIndex;
            break;
        end
    end
    if ~foundMatch
        distinctLabels = [distinctLabels currentLabels];
        labelCount = labelCount + 1;
        labelIndices(trialIndex) = labelCount;
    end
end