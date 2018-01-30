function [repairedData, coefficients, grossErrors] = repair_corrupted_data(dictionary, corruptedData)
% repair_missing_data.m
%
%   Given a set of complete data vectors, and another set of vectors with
%   missing entries, use L1-minimization to repair the incomplete vectors.
%   This function required the CVX package for semidefinite programming.
%
%
% Inputs:
%   dictionary           - a matrix whose columns are data vectors that are
%                          used as overcomplete basis to represent other
%                          vectors. note that some or all of these columns
%                          may contain gross errors.
%   corruptedData        - a matrix whose columns are data vectors with
%                          some elements that may have gross errors. these
%                          vectors will be repaired using the dictionary.
%
% Outputs:
%   repairedData        - a matrix whose columns are repaired versions of
%                         the given incomplete data vectors
%   coefficients        - a matrix whose ith column is  a list of
%                         coefficients used to represent the ith vector in
%                         in corrupted data as a linear combination of
%                         vectors in the dictionary
%   grossErrors         - a matrix whose ij-th entry is true if the ij-th
%                         entry of corruptedData was corrupted by gross
%                         errors.
% Dependencies:
%   CVX package
%
% Nov. '07  Shankar Rao -- srrao@uiuc.edu

% Copyright 2007, University of Illinois. All rights reserved.

VERBOSE = false;
GROSS_ERROR_THRESHOLD = 0.5;

[dimensionCount, sampleCount] = size(dictionary);
corruptedSampleCount = size(corruptedData,2);


repairedData = corruptedData;
grossErrors = false(size(corruptedData));
coefficients = zeros(sampleCount+dimensionCount, corruptedSampleCount);
normalizedData = dictionary ./ repmat(sqrt(sum(dictionary.^2, 1)), dimensionCount, 1);
I = eye(dimensionCount);
X = [ normalizedData I];
for sampleIndex = 1:corruptedSampleCount
    y = corruptedData(:, sampleIndex);

    if VERBOSE,
        disp(sprintf('Repairing sample %d of %d', sampleIndex, corruptedSampleCount));
    end

    state = cvx_quiet(true);
    cvx_begin
        variable c(sampleCount+dimensionCount);
        minimize(norm(c,1));
            subject to
                X * c == y;
    cvx_end
    cvx_quiet(state);

    yhat = X(:, 1:sampleCount)*c(1:sampleCount);
    coefficients(:, sampleIndex) = c;
    errors = c(sampleCount+1:end);
    grossErrors(:, sampleIndex) = abs(errors) > GROSS_ERROR_THRESHOLD;
    repairedData(grossErrors(:, sampleIndex), sampleIndex) = yhat(grossErrors(:, sampleIndex));
end
