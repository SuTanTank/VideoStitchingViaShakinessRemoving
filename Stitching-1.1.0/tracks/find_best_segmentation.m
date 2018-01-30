function [bestLabels, bestEpsilon] = find_best_segmentation(result, sequenceData,  groupCount, epsilon)

% find_best_segmentation.m
%
%   Given a known group count, choose the segmentation in a list of segmentations
%   that best codes the data. The approach uses a modified voting scheme to
%   pick a segmentation that has the correct number of groups. If no such
%   segmentation exists, than the segmentation with the minimum penalized
%   coding length is returned, along with a warning.
%
% Inputs:
%   result              - a structure containing the results of applying
%                           greedy algorithm to a particular motion
%                           sequence for a range of epsilons.
%   sequenceData        - a matrix where each column is a vector obtained
%                           from a tracked feature in the motion sequence.
%   groupCount          - the known number of groups.
%   epsilon             - a list containing the range of epsilons for which these results
%                           were obtained.
%
% Outputs:
%   bestLabels          - a labeling of the motion data that is the "best"
%                           segmentation.
%   bestEpsilon         - a choice of epsilon that will cause ALC to
%                         produce this "best" segmentation.
% Dependencies:
%   distinct_labels, total_coding_length
%
% Mar. '08  Shankar Rao -- srrao@uiuc.edu
%
% Copyright 2008, University of Illinois. All rights reserved.

OUTLIER_GROUP_COUNT = 4;

[sampleCount, trialCount] = size(result.labels);
properGroupings = [];
isSingleGrouping = false(1,trialCount);

for trialIndex = 1:trialCount
    currentLabels = result.labels(:, trialIndex);
    distribution = histc(currentLabels, 1:max(currentLabels));
    if (sum(distribution > OUTLIER_GROUP_COUNT) == groupCount)
        properGroupings = [properGroupings trialIndex];
        % Label all outliers as one group, relabel inliers
        % modified 10-24-07
        outlierLabels = find(distribution <= OUTLIER_GROUP_COUNT);
        inlierLabels = find(distribution > OUTLIER_GROUP_COUNT);
        outlierIndices = find(ismember(currentLabels, outlierLabels));
        currentLabels(outlierIndices) = -1;
        for labelIndex = 1:length(inlierLabels)
            currentLabels(result.labels(:, trialIndex) == inlierLabels(labelIndex)) = labelIndex;
        end
        currentLabels(outlierIndices) = groupCount+1;
        result.labels(:, trialIndex) = currentLabels;
    elseif (max(currentLabels) == 1)
        isSingleGrouping(trialIndex) = true;
    end
end


if isempty(properGroupings)
    warning('No segmentation has the correct number of groups. Returning the one with minimal penalized coding length.');
    [ignore winnerIndex] = min(result.penalty);
    bestLabels = result.labels(:, winnerIndex);
    bestEpsilon = epsilon(winnerIndex);
else
    % ignore any epsilon that results in all of the samples being grouped
    % together.
    nonSingleGroupings = find(~isSingleGrouping);

    % only epsilons that result in the proper number of group can be
    % considered as candidate groupings.
    properGroupingsCount = length(properGroupings);

    [distinctLabels labelIndices] = distinct_labels(result.labels(:, properGroupings));
    candidateCount = size(distinctLabels, 2);
    codingLengths = zeros(candidateCount, length(nonSingleGroupings));

    % however, all non single-group segmentations get a vote
    for trialIndex = 1:length(nonSingleGroupings)
        for candidateIndex = 1:candidateCount
            codingLengths(candidateIndex, trialIndex) = ...
                total_coding_length(sequenceData, distinctLabels(:, candidateIndex), epsilon(nonSingleGroupings(trialIndex)), false);
        end
    end

    [ignore votes] = min(codingLengths, [] , 1);
    voteCounts = histc(votes, 1:candidateCount);
    [ignore winnerIndex] = max(voteCounts);
    bestLabels = distinctLabels(:, winnerIndex);
    winnerLabelIndices = find(labelIndices == winnerIndex);
    bestEpsilon = epsilon(uint32(floor(median(winnerLabelIndices))));
end

