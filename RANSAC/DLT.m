% Fit homography
function H = DLT(varargin)

    % Check the input parameters
    [x1, x2] = checkargs(varargin(:));
    
    % Normalise each set of points 
    [x1, T1] = NormalisePoints(x1);
    [x2, T2] = NormalisePoints(x2);
    
    Npts = length(x1);
    A = zeros(3*Npts,9);
    
    O = [0 0 0];
    for n = 1:Npts
        X = x1(:,n)';
        x = x2(1,n); 
        y = x2(2,n); 
        w = x2(3,n);
        A(3*n-2,:) = [  O  -w*X  y*X];
        A(3*n-1,:) = [ w*X   O  -x*X];
        A(3*n  ,:) = [-y*X  x*X   O ];
    end
    
    % SVD decompose
    [U,D,V] = svd(A,0); 
    
    % Extract homography
    H = reshape(V(:,9),3,3)';
    
    % Denormalise
    H = T2\H*T1;
    

% Check the input parameters
function [x1, x2] = checkargs(arg)
    
    if length(arg) == 2
        x1 = arg{1};
        x2 = arg{2};
        if ~all(size(x1)==size(x2))
            error('x1 and x2 must have the same size');
        elseif size(x1,1) ~= 3
            error('x1 and x2 must be 3xN');
        end

    elseif length(arg) == 1
        if size(arg{1},1) ~= 6
            error('Single argument x must be 6xN');
        else
            x1 = arg{1}(1:3,:);
            x2 = arg{1}(4:6,:);
        end
    else
        error('Wrong number of arguments supplied');
    end
    
    
    