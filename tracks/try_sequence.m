function result = try_sequence(sequenceName, sequenceData, epsilon, trueLabels, labels)

% try_sequence
%
%   Try segmenting a set of motion data using the greedy algorithm over a
%   range of epsilons.
% 
% Inputs:
%   sequenceName - the name of this motion sequence to be tested                 
%   sequenceData - a matrix where each column is a processed feature vector
%                   of the motion sequence.
%   epsilon      - the error tolerance
%   trueLabels   - the apriori known labelling of the data. this is used to
%    (optional)     compute the misclassification error for each trial.
%   labels       - labels that we computed before. if this is passed to the
%    (optional)     function, then we are just trying to reconstruct the result structure
%
% Outputs:
%   result - a structure containing the results of applying
%             greedy algorithm to a particular motion
%             sequence for a range of epsilons.
% Dependencies:
%   coding_seg, compare_labels, total_coding_length
%
% Mar. '08  Shankar Rao -- srrao@uiuc.edu
%
% Copyright 2008, University of Illinois. All rights reserved.

% tic
VERBOSE = true;
trialCount = length(epsilon);

result.name = sequenceName;
[dimensionCount, sampleCount] = size(sequenceData);
if nargin < 5
    RECOMPUTE_LABELS = true;
    result.labels = zeros(sampleCount, trialCount);
else
    RECOMPUTE_LABELS = false;
    result.labels = labels;
end

COMPUTE_ERROR = (nargin >= 4);

result.groupCounts = zeros(1, trialCount);
result.coding = zeros(1, trialCount);
result.penalty = zeros(1, trialCount);
result.error = ones(1, trialCount);

for trialIndex=1:trialCount
    if RECOMPUTE_LABELS,
       [result.labels(:, trialIndex), codingLengths] = coding_seg(sequenceData, epsilon(trialIndex), false);
    end
    computedCodingLength = total_coding_length(sequenceData, result.labels(:, trialIndex)', epsilon(trialIndex), false);
    oneGroupCodingLength = total_coding_length(sequenceData, ones(1,sampleCount), epsilon(trialIndex), false);
    if (computedCodingLength > oneGroupCodingLength)
        if VERBOSE,
            disp(sprintf('For epsilon=%f computed coding length=%f > coding length of data as one group=%f', epsilon(trialIndex), computedCodingLength, oneGroupCodingLength));
        end
        result.labels(:, trialIndex) = ones(sampleCount, 1);
        result.coding(trialIndex) = oneGroupCodingLength;
    else
        result.coding(trialIndex) = computedCodingLength;
    end
    result.groupCounts(trialIndex) = max(result.labels(:, trialIndex));
    if COMPUTE_ERROR && (result.groupCounts(trialIndex) <= 50)
        result.error(trialIndex) = compare_labels(trueLabels, result.labels(:,trialIndex)');
    end
    result.penalty(trialIndex) = result.coding(trialIndex) + (dimensionCount*sampleCount*log2(epsilon(trialIndex)));
    if VERBOSE
        disp(sprintf('  --> Trial %d, # of motions=%d, epsilon=%f, error=%f, coding length=%f, pen. coding length=%f', ...
            trialIndex, result.groupCounts(trialIndex),epsilon(trialIndex), result.error(trialIndex), result.coding(trialIndex), result.penalty(trialIndex)));
    end
end

% result.time = toc;
% if VERBOSE
%     disp(sprintf('     Time Elapsed: %d seconds', result.time));
% end

