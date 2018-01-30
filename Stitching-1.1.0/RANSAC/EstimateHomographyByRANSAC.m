
% Estimate 2D homography using RANSAC
function [H, inliers] = EstimateHomographyByRANSAC(x1, x2, t)

    if ~all(size(x1)==size(x2))
        error('Data sets x1 and x2 must have the same dimension');
    end
    
    [rows,npts] = size(x1);
    if rows~=2 && rows~=3
        error('x1 and x2 must have 2 or 3 rows');
    end
    
    if npts < 4
        error('Must have at least 4 points to fit homography');
    end
    
    if rows == 2    
        x1 = [x1; ones(1,npts)];
        x2 = [x2; ones(1,npts)];        
    end
    
%     % compute centering of points % by Tan
%     if videoSize(1) == 0
%         scoreX1 = 0;
%         scoreX2 = 0;
%     else
%         mid = videoSize / 2;
%         p1 = x1(1:2, :);
%         p2 = x2(1:2, :);
%         p1x = p1(1, :);
%         p1y = p1(2, :);
%         p2x = p2(1, :);
%         p2y = p2(2, :);
%         p1x = abs(p1x - mid(1))/mid(1);
%         p1y = abs(p1y - mid(2))/mid(2);
%         p2x = abs(p2x - mid(1))/mid(1);
%         p2y = abs(p2y - mid(2))/mid(2);
%         scoreX1 = max([p1x;p1y]);
%         scoreX2 = max([p2x;p2y]);    
%     end
%     
    
    % Normalise each set of points
    [x1, T1] = NormalisePoints(x1);
    [x2, T2] = NormalisePoints(x2);
    
    % Minimum No of points needed to fit a homography.
    s = 4;  
    
    DLTfn = @DLT;
    errorfn    = @homogError;
    degenfn   = @checkDegenerate;
    
    [H, inliers] = RANSAC([x1; x2], DLTfn, errorfn, degenfn, s, t);
    H = DLT(x1(:,inliers), x2(:,inliers));
    
    % Denormalise
    H = T2\H*T1;    
end


% reprojection error.
function [inliers, score, H] = homogError(H, x, t)
    
    x1 = x(1:3,:);
    x2 = x(4:6,:);    
    
    % Calculate, in both directions, the transfered points    
    Hx1    = H*x1;
    invHx2 = H\x2;
    
    % Normalise, ensure the scalar is 1
    x1     = hnormalise(x1);
    x2     = hnormalise(x2);     
    Hx1    = hnormalise(Hx1);
    invHx2 = hnormalise(invHx2); 
    
    d2 = sum((x1-invHx2).^2)  + sum((x2-Hx1).^2);
    inliers = find(abs(d2) < t);    
    
    % Su Tan
    % score = min(max(x1(1, inliers)) - min(x1(1, inliers)) + max(x1(2, inliers)) - min(x1(2, inliers)), ...
        %max(x2(1, inliers)) - min(x2(1, inliers)) + max(x2(2, inliers)) - min(x2(2, inliers)));
    score = 1;
    score = score * length(inliers);
end   
    
% normalise 
function nx = hnormalise(x)
    
    [rows,npts] = size(x);
    nx = x;

    % Find the indices of the points that are not at infinity
    finiteind = find(abs(x(rows,:)) > eps);

    if length(finiteind) ~= npts
        disp('Some points are at infinity');
    end

    % Normalise points not at infinity
    for r = 1:rows-1
        nx(r,finiteind) = x(r,finiteind)./x(rows,finiteind);
    end
        nx(rows,finiteind) = 1;
end
    
% check if the three points are colinear
function r = iscolinear(p1, p2, p3, flag)

    if nargin == 3   % Assume inhomogeneous coords
	flag = 'inhomog';
    end
    
    if ~all(size(p1)==size(p2)) || ~all(size(p1)==size(p3)) || ...
        ~(length(p1)==2 || length(p1)==3)                              
        error('points must have the same dimension of 2 or 3');
    end
    
    if length(p1) == 2    
        p1(3) = 1; p2(3) = 1; p3(3) = 1;
    end

    if flag(1) == 'h'
        r =  abs(dot(cross(p1, p2),p3)) < eps;
    else
        r =  norm(cross(p2-p1, p3-p1)) < eps;
    end
    
end    
% Check whether any 3 of the 4 points in each set is colinear. 
function r = checkDegenerate(x)

    x1 = x(1:3,:);
    x2 = x(4:6,:);    
    
    r = iscolinear(x1(:,1),x1(:,2),x1(:,3)) | iscolinear(x1(:,1),x1(:,2),x1(:,4)) | ...
    iscolinear(x1(:,1),x1(:,3),x1(:,4)) | iscolinear(x1(:,2),x1(:,3),x1(:,4)) | ...
    iscolinear(x2(:,1),x2(:,2),x2(:,3)) | iscolinear(x2(:,1),x2(:,2),x2(:,4)) | ...
    iscolinear(x2(:,1),x2(:,3),x2(:,4)) | iscolinear(x2(:,2),x2(:,3),x2(:,4));
end
