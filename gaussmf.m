function y = gaussmf(x, params)
%GAUSSMF Gaussian curve membership function.
%   GAUSSMF(X, PARAMS) returns a matrix which is the Gaussian
%   membership function evaluated at X. PARAMS is a 2-element vector
%   that determines the shape and position of this membership function.
%   Specifically, the formula for this membership function is:
%
%   GAUSSMF(X, [SIGMA, C]) = EXP(-(X - C).^2/(2*SIGMA^2));
%   
%   For example:
%
%       x = (0:0.1:10)';
%       y1 = gaussmf(x, [0.5 5]);
%       y2 = gaussmf(x, [1 5]);
%       y3 = gaussmf(x, [2 5]);
%       y4 = gaussmf(x, [3 5]);
%       subplot(211); plot(x, [y1 y2 y3 y4]);
%       y1 = gaussmf(x, [1 2]);
%       y2 = gaussmf(x, [1 4]);
%       y3 = gaussmf(x, [1 6]);
%       y4 = gaussmf(x, [1 8]);
%       subplot(212); plot(x, [y1 y2 y3 y4]);
%       set(gcf, 'name', 'gaussmf', 'numbertitle', 'off');
%
%   See also DSIGMF,EVALMF, GAUSS2MF, GBELLMF, MF2MF, PIMF, PSIGMF, SIGMF,
%   SMF, TRAPMF, TRIMF, ZMF.

%   Roger Jang, 6-29-93, 10-5-93.
%   Copyright (c) 1994-98 by The MathWorks, Inc.
%   $Revision: 1.14 $  $Date: 1998/08/06 13:57:19 $

if nargin ~= 2
    error('Two arguments are required by the Gaussian MF.');
elseif length(params) < 2
    error('The Gaussian MF needs at least two parameters.');
elseif params(1) == 0,
    error('The Gaussian MF needs a non-zero sigma.');
end

sigma = params(1); c = params(2);
y = exp(-(x - c).^2/(2*sigma^2));