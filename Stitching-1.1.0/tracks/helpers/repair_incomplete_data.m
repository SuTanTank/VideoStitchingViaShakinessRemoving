function [repairedData, coefficients] = repair_incomplete_data(completeData, incompleteData, missingEntries, method)
% repair_incomplete_data.m
%
%   Given a set of complete data vectors, and another set of vectors with
%   missing entries, use L1-minimization to repair the incomplete vectors.
%   This function requires the CVX package 
%
% Inputs:
%   completeData        - a matrix whose columns are complete data vectors
%   incompleteData      - a matrix whose columns are incomplete data
%                         vectors. The missing entries can be set to
%                         anything, because they will be ignored
%   missingEntries      - A matrix whose ij-th entry is 1 if the ith
%                         element of the jth incomplete data vector is missing, 
%                         and 0 otherwise
%   method (optional)   - if the string 'completion'(default), perform
%                           L1-minimization by removing the rows with
%                           missing entries. if the string 'correction',
%                           perform L1-minimization by adding columns of
%                           the identity matrix corresponding to the
%                           missing entries to the complete data matrix and
%                           perform error correction.
%
% Outputs:
%   repairedData        - a matrix whose columns are repaired versions of
%                         the given incomplete data vectors
%   coefficients        - A matrix whose ith column is the coefficients
%                         required to express the ith incomplete data
%                         vector as a linear combination of the vectors in 
%                         completeData  
% Dependencies:
%   CVX package
%
% Nov. '07  Shankar Rao -- srrao@uiuc.edu
%
% Copyright 2007, University of Illinois. All rights reserved.

VERBOSE = false;

if (nargin == 3)
    method = 'completion';
end

[dimensionCount, sampleCount] = size(completeData);
incompleteDataCount = size(incompleteData, 2);


X = completeData ;
repairedData = incompleteData;
coefficients = nan(sampleCount, incompleteDataCount);

for sampleIndex = 1:incompleteDataCount

    y = incompleteData(:, sampleIndex);

    entries = find(missingEntries(:, sampleIndex));
    entryCount = length(entries);
    if VERBOSE,
        disp(sprintf('Repairing sample %d of %d with %d missing entries.', sampleIndex, incompleteDataCount, entryCount));
    end
    if strcmp(method, 'completion')
        existingEntries = setdiff(1:dimensionCount, entries);

        % Set up linear program variables
        X_tilde = X(existingEntries, :);
        y_tilde = y(existingEntries);

        if (size(X_tilde, 1) > size(X_tilde,2))
            disp('Warning: Not enough complete data for reconstruction');
            repairedData(entries, sampleIndex) = 0;
            continue;
        end
        % Solve linear program using cvx package
        state = cvx_quiet(true);
        cvx_begin
          variable c(sampleCount);
          minimize(norm(c,1));
          subject to
            X_tilde * c == y_tilde;
        cvx_end
        cvx_quiet(state);

        yhat = X*c;        
        yhat(existingEntries) = y_tilde;
%    elseif strcmp(method, 'correction')
%         X = X ./ repmat(sqrt(sum(completeData.^2, 1)), dimensionCount,1);
%         I = eye(dimensionCount);
%         phiI = I(:, entries);
%         X_tilde = [X phiI];
%         y_tilde = y;
%         y_tilde(entries) = 0;
% 
% 
%         state = cvx_quiet(true);
%         cvx_begin
%             variable c(sampleCount+missingEntryCount);
%             minimize(norm(c,1));
%             subject to
%                 X_tilde * c == y_tilde;
%         cvx_end
%         cvx_quiet(state);
%         
%         yhat = X_tilde*c - phiI*c(end-entryCount+1:end);
    end
    repairedData(entries, sampleIndex) = yhat(entries);
    coefficients(:, sampleIndex) = c;
end
