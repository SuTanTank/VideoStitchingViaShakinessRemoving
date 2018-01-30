function [sampleLabels, groupCodingLengths] = coding_seg(W, epsilon, affine )

% coding_seg.m
%
% Implements an agglomerative segmentation algorithm which attempts to
% minimize the number of bits needed to code the segmented vectors upto
% distortion epsilon^2. This method is described in detail in Ma, et. al.,
% Segmentation of Multivariate Mixed Data via Lossy Coding and Compression,
% PAMI 2007.
%
% Run test_coding_seg.m, included in this package, for a demonstration. 
%
% Inputs:
%   W       - [ w1 | w2 | ... | wm ], the data matrix. 
%               The columns are the input data vectors. 
%
%   epsilon - sqrt of the allowable distortion, 
%               E[ \| \hat{w_i} - w_i \|^2 ]
%
%   affine  - boolean, if true the algorithm uses the coding length for
%               nonzero-mean data, derived in Appendix II of the paper.
%
% Outputs:
%   sampleLabels - 1 x m list whose i-th entry is the group assignment of
%                   the i-th data vector w_i. Groups are indexed
%                   sequentially, starting from 1. 
%
% Dependencies:
%   coding_length.m
%
% Feb '06, last revision Jan '07
%   Questions? John Wright -- jnwright@uiuc.edu

% Copyright 2007, University of Illinois. All rights reserved. 

    VERBOSE = 0;

    epsilon2 = epsilon^2;

    % n - dimension 
    % m - number of samples
    [n,m] = size(W);

    % initially treat each sample as its own group
    sampleLabels = 1:m;
    curGroupCount = m;

    % compute the coding length of each of the initial groups
    squaredNormW = sum(W.^2,1);
    if affine
        groupCodingLengths = n * log2(1 + (squaredNormW/epsilon2))/ 2 + log2(m);
    else
        groupCodingLengths = (n + 1) * log2(1 + (n / epsilon2)*squaredNormW)/ 2 + log2(m);
    end

    if VERBOSE
        disp('   Computing initial table.');
    end

    % entropyChange is a matrix whose ij-th element is the entropy lost by
    % merging groups i and j into one group.
    delta_L_matrix = ones(m);

    for groupIndex1 = 1 : m-1
        for groupIndex2 = groupIndex1+1 : m
            L1 = groupCodingLengths(groupIndex1);
            L2 = groupCodingLengths(groupIndex2);
            L_union = coding_length(W(:, [groupIndex1 groupIndex2]),epsilon2,m,affine);
            delta_L_matrix(groupIndex1, groupIndex2) = L_union - L1 - L2;
        end;
    end;

    if VERBOSE
        disp('   Starting merging process');
    end

    while curGroupCount > 1 

        if mod( curGroupCount, 50 ) == 0 && VERBOSE
            disp(['   Group count: ' num2str(curGroupCount)]);
        end

        % Find the best two groups to merge together
        [temp minIndex]=min(delta_L_matrix);
        [min_delta_L minIndex2]=min(temp);
        minIndex1=minIndex(minIndex2);

        % If this does not decrease the coding length, we are done
        if ( min_delta_L > 0 )
            break;
        end;

        % Merging groups with labels minIndex1 and minIndex2
        if (minIndex1 > minIndex2)
            temp = minIndex1;
            minIndex1 = minIndex2;
            minIndex2 = temp;
        end

        if (minIndex2 ~= curGroupCount)
            % First make the group to merge have the last label

            % Change the sample labels of the group minIndex2 to be the last
            % label.
            sampleLabels(sampleLabels == curGroupCount) = -1;
            sampleLabels(sampleLabels == minIndex2) = curGroupCount;
            sampleLabels(sampleLabels == -1) = minIndex2;

            % Swap the entropy of groups minIndex2 and the last group.        
            [groupCodingLengths(minIndex2), groupCodingLengths(curGroupCount)] = ...
                swap(groupCodingLengths(minIndex2), groupCodingLengths(curGroupCount));

            delta_L_matrix = swap_matrix(delta_L_matrix, minIndex2, curGroupCount);

        end

        % Now we can assume the last group is the group to be merged    
        sampleLabels(sampleLabels == curGroupCount) = minIndex1;
        groupCodingLengths(curGroupCount) = [];
        groupCodingLengths(minIndex1) = coding_length(W(:, sampleLabels == minIndex1),epsilon2,m,affine);
        delta_L_matrix = delta_L_matrix(1:end-1, 1:end-1);

        curGroupCount = curGroupCount - 1;    

        for groupIndex2 = [1:minIndex1-1 minIndex1+1:curGroupCount],
            i = min(minIndex1, groupIndex2); j = max(minIndex1, groupIndex2); 
                L1 = groupCodingLengths(i);
                L2 = groupCodingLengths(j);
                L_union = coding_length(W(:, (sampleLabels == i) | (sampleLabels == j)) ,epsilon2,m,affine);
                delta_L_matrix(i, j) = L_union - L1 - L2;            
        end;
    end;

    %codingLength = sum(groupCodingLengths);

    if VERBOSE
        disp(['   Final group count: ' num2str(curGroupCount)]);
    end

end
% If A is an upper triangular matrix where A(m,n) = f(x_m,x_n) for m < n
% relabel so that the roles of x_i and x_j are swapped.
function A = swap_matrix(A, i, j)
[A(1:i-1, i), A(1:i-1, j)] = swap(A(1:i-1, i), A(1:i-1, j));
[A(i, i+1:j-1), A(i+1:j-1, j)] = swap(A(i, i+1:j-1), A(i+1:j-1, j));
[A(i, j+1:end), A(j, j+1:end)] = swap(A(i, j+1:end), A(j, j+1:end));
end

% swap the values of x and y
function [y, x] = swap(x,y)
end