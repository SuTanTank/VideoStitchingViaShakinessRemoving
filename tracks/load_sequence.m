function [data trueLabels mask] = load_sequence(sequenceName, normalizeData)

% load_sequence.m
%
%   Load a particular motion sequence.
%
% Inputs:
%   sequenceName                    - The name of the motion sequence to
%                                       load.
%
% Outputs:
%   trueLabels                      - A list of labels for the correct segmentation of the data
%   data                            - A matrix of the raw tracked feature correspondences
%   normalizedData                  - A matrix of the tracked feature
%                                      correspondences normalized to range between -1 and 1.
%   mask                            - contains the 'mask' for each point in each frame. More precisely:
%                                     mask(p,f)==1: point p is visible in the frame f
%                                     mask(p,f)==0: point p is missing in
%                                     the frame  
% Dependencies:
%   none
%
% Mar. '08  Shankar Rao -- srrao@uiuc.edu

% Copyright 2008, University of Illinois. All rights reserved.

if nargin < 2
    normalizeData = false;
end

SEQUENCE_FOLDER_NAME = 'datasets';
fileName = [ SEQUENCE_FOLDER_NAME '/' sequenceName '.mat'];
load(fileName);

trueLabels = s;
if normalizeData
    data = x;
else 
    data = y;
end

if nargout >= 3
    if exist('m', 'var') && nargout >= 3
        mask = m;
    else
        error(sprintf('No mask found for sequence %s', sequenceName));
    end
end