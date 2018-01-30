function [processedData complete isGrossError ] = process_sequence(rawData, projection, mode, mask)

% process_sequence
%
%       process tracked feature correspondences. the raw data is assumed to
%       be in a 2xFxP matrix, where F is the number of frames and P is the
%       number of points. For each tracked feature, it coordinates in all
%       frames are stacked together to obtain a 2Fx1 vector. additionally, 
%       projection onto a low-dimensional subspace, or dealing with incomplete or
%       corrupted trajectories may also be required.
%
% Inputs:
%   rawData             - A matrix of the raw tracked feature
%                           correspondences.
%
%   projection          - if true, project the feature vector onto a
%                           5-dimensional space using PCA. if false, omit
%                           any projection. If equal to the string
%                           "sparse", project the feature vectors ont a
%                           d-dimension space using PCA, where d is the
%                           smallest integer that satisfies d >
%                           2*k*lg(N/d).
%
%   mode                - if equal to 'incomplete', then fill in incomplete trajectories 
%    (optional)           specified by mask. if equal to 'corrupted', detected and repair 
%                         corrupted trajectories. otherwise, the mode should be missing 
%                         or equal to 'clean'. 
%   mask                - if mode equals 'incomplete', then this is a required PxF matrix 
%    (optional)           whose ijth entry is 1 if the ith point is visible in the jth frame, 
%                         and 0 otherwise. this matrix is ignored for all other modes.        
%                  
% Outputs:
%   processedData        - the processed motion data.
%
%   isGrossError         - if mode equals 'corrupted', then this is a 2FxP
%                          matrix whose ij-th entry is true if the ith
%                          entry of the jth trajectory was corrupted by a
%                          gross error.
% Dependencies:
%   repair_incomplete_data, repair_corrupted_data
%
% Mar. '08  Shankar Rao -- srrao@uiuc.edu

% Copyright 2008, University of Illinois. All rights reserved.


VERBOSE = false;

[ignore sampleCount frameCount] = size(rawData);

data = zeros(2*frameCount, sampleCount);
if nargin < 3
    mode = 'clean';
end

if strcmp(mode, 'incomplete')
    isComplete = all(mask,2);
    completeSampleCount = sum(isComplete);
    incompleteSampleCount = sum(~isComplete);

    completeData = zeros(2*frameCount, completeSampleCount);
    incompleteData = zeros(2*frameCount, incompleteSampleCount);
    incompleteMask = ones(2*frameCount, incompleteSampleCount);
    cSampleIndex = 1;     % complete sample index
    icSampleIndex = 1;    % incomplete sample index
    for sampleIndex=1:sampleCount
        if isComplete(sampleIndex)
            completeData(:, cSampleIndex) = reshape(rawData(1:2, sampleIndex, :), 2*frameCount,1);
            cSampleIndex = cSampleIndex + 1;
        else
            incompleteData(:, icSampleIndex) = reshape(rawData(1:2, sampleIndex, :), 2*frameCount,1);
            incompleteMask(1:2:2*frameCount-1, icSampleIndex) = mask(sampleIndex, :)';   % ??2F * P???
            incompleteMask(2:2:2*frameCount, icSampleIndex) = mask(sampleIndex, :)';
            icSampleIndex = icSampleIndex + 1;
            if VERBOSE,
              disp(sprintf('Completed sample %d of %d with %d missing points', icSampleIndex, incompleteSampleCount, sum(~mask(sampleIndex,:))));
            end
        end

    end
        
    repairedData = repair_incomplete_data(completeData, incompleteData, ~incompleteMask);
    data(:, isComplete) = completeData;
    data(:, ~isComplete) = repairedData;
    complete = data;
elseif strcmp(mode, 'clean') || strcmp(mode, 'corrupted')
    for sampleIndex=1:sampleCount
        data(:, sampleIndex) = reshape(rawData(1:2, sampleIndex, :), 2*frameCount,1);
    end
else
    error('Unknown mode for processing motion sequence data.');
end

if strcmp(mode, 'corrupted')
    corruptedData = data;
    isGrossError = false(size(data));
    for sampleIndex = 1:sampleCount
        [data(:, sampleIndex), ignore, isGrossError(:, sampleIndex)] = ...
            repair_corrupted_data(corruptedData(:, [1:sampleIndex-1 sampleIndex+1:sampleCount]), corruptedData(:, sampleIndex));
        if VERBOSE,
            disp(sprintf('Repaired sample %d of %d, %d errors found', sampleIndex, sampleCount, sum(isGrossError(:, sampleIndex))));
        end
    end
end


% projecting data into lower dimensional space
if projection,
    if strcmp(projection, 'sparse')
        dimensionCount = 1;
        while (dimensionCount < 8*log2(2*frameCount/dimensionCount))
            dimensionCount = dimensionCount + 1;
        end
        if VERBOSE
            disp(sprintf('Calculated dimension = %d', dimensionCount));
        end
    else
        dimensionCount = 5;
    end
    [U S] = svd(data, 0);
    basis = U(:, 1:dimensionCount);
    processedData = basis'*data;
else
    dimensionCount = size(data, 1);
    processedData = data;
end