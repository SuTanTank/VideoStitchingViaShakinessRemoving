function [bestMapping, errorProbability] = relabel_samples(aprioriSampleLabels, computedSampleLabels, subspaceDimensions)
%function [bestMapping, errorProbability, groupError] = relabel_samples(aprioriSampleLabels, computedSampleLabels, subspaceDimensions)
%
% Attempts to find a good mapping between two sets of sample labels based on their
% grouping alone. A valid permutation vector will be returned in all cases.
% 
% The function is restricted to only permute subspaces of equal dimensions.
%
% Sep. '05  Allen Yang -- yang@eecs.berkeley.edu, Shankar Rao -- srrao@uiuc.edu

% Copyright 2005, University of Illinois. All rights reserved.

if sort(subspaceDimensions,'descend')~=subspaceDimensions
    error('The parameter subspaceDimensions should be in descending order');
end

% Eliminate outlier labels
inlierIndex=find((aprioriSampleLabels~=-1).*(computedSampleLabels~=-1));
aprioriSampleLabels = aprioriSampleLabels(inlierIndex);
computedSampleLabels = computedSampleLabels(inlierIndex);

groupCount = length(subspaceDimensions);
computedSampleFrequencies = histc(computedSampleLabels, 1:groupCount);
sampleCount = length(aprioriSampleLabels);
errorProbability = 1;
newSampleLabels = zeros(1,sampleCount);

% Generate all legal permutations. No cross permutation of different
% dimensions
allPermutations = [];
for dimensionIndex=max(subspaceDimensions):-1:1
    indices = find(subspaceDimensions==dimensionIndex);
    if length(indices)>0
        if isempty(allPermutations)
           allPermutations = perms(indices); 
        else
            previousPermutation = allPermutations;
            previousPermutationSize = size(previousPermutation,1);
            partialPermutation = perms(indices);
            allPermutations = [previousPermutation,repmat(partialPermutation(1,:),[previousPermutationSize,1])];
            for permutationIndex=2:size(partialPermutation,1)
                allPermutations = [allPermutations; ...
                    previousPermutation,repmat(partialPermutation(permutationIndex,:),[previousPermutationSize,1])];
            end
        end
    end
end


embedded = groupCount*(computedSampleLabels - 1) + (aprioriSampleLabels - 1);
histogram = histc(embedded, 0:groupCount^2 - 1);
for permutationIndex = 1 : size(allPermutations, 1)
    correct = 0;
    mapping = allPermutations(permutationIndex, :);

    correct = sum(histogram(groupCount*(0:groupCount-1) + mapping)) / sampleCount;
    
    if (correct > 1 - errorProbability)
        errorProbability = 1 - correct;
        bestMapping = mapping;
    end
end