function [bestError bestMap] = compare_labels(correctLabels, computedLabels)

% compare_labels
%
%   Given two labelings of a set of data, find a map between the
%   labelings that minimizes the number of samples that they assign
%   different labels. This version assumes the correct labeling has a small
%   number of groups (<= 5), and does an exhaustive search to find
%   the best map. If the correct number of groups is larger, the result
%   from relabel_samples (which will likely be incorrect) is run.
%
% Inputs:
%   correctLabels - the true labels of some data.
%   computedLabels - the labels computed by some algorithm for that data.
%
% Outputs:
%   bestError - the misclassification rate obtained by this relabelling.
%   bestMap - the map that converts correct labels to computed labels
%   
% Sep. '07  Shankar Rao -- srrao@uiuc.edu

% Copyright 2007, University of Illinois. All rights reserved.

correctLabels = correctLabels(:)';
computedLabels = computedLabels(:)';

groupCount = max(correctLabels);
computedGroupCount = max(computedLabels);
sampleCount = length(correctLabels);

maxGroupCount = max(groupCount, computedGroupCount);

if groupCount <= 5
    bestError = 1;
    bestMap = [1:groupCount];
    candidateLabels = nchoosek(1:maxGroupCount, groupCount);
    permCount = factorial(groupCount);
    for candidateIndex = 1:size(candidateLabels,1)
        maps = perms(candidateLabels(candidateIndex,:));
        for permIndex = 1:permCount
            map = maps(permIndex,:);
            currentError = mean(map(correctLabels) ~= computedLabels);
            if currentError < bestError
                bestError = currentError;
                bestMap = map;
            end
        end
    end
else
    [bestMap, bestError] = relabel_samples(correctLabels, computedLabels, ones(1,groupCount));
end