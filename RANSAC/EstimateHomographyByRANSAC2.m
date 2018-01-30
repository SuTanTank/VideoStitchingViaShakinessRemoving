% Estimate 2D homography using RANSAC
function [H, inliers, score] = EstimateHomographyByRANSAC2(x1, x2, t, index, validMask)

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
    
    do = 1;
    
    while do
    
        [~, inliers] = RANSAC2([x1; x2], DLTfn, errorfn, degenfn, s, t, index, validMask);

    %     len = length(inliers);
        inliers_less = inliers;
    %     for k = 1:len
    %         if validMask(index(inliers(k))) == 0
    %             inliers_less(k) = 0; 
    %         end
    %     end
    %     inliers_less(inliers_less == 0) = [];

        H = DLT(x1(:,inliers_less), x2(:,inliers_less));
        [inliers, score, H] = homogError(H, [x1; x2], t, index, validMask);
        % Denormalise
        H = T2\H*T1;    

        x_index = floor((index(inliers) - 1) / 10) + 1;
        y_index = index(inliers) - x_index * 10;
        dis = (max(x_index) - min(x_index) + 1) * (max(y_index) - min(y_index) + 1) / 100;

        voting = histcounts([index(inliers); 1; 100], 100);  
        voting(1) = voting(1) - 1;
        voting(100) = voting(100) - 1;

        voting_all = histcounts([index; 1; 100], 100);  
        voting_all(1) = voting_all(1) - 1;
        voting_all(100) = voting_all(100) - 1;
        % Tan Su
        voting(voting < voting_all * 0.5) = 0;
        voting(voting > 0) = 1;
    %     u = sum(validMask' | voting);
    %     i = sum(validMask' & voting);
    %     fprintf('%d/%d\n', i, u);

        score = sum(voting & validMask') / sum(voting | validMask');
        
        if score < 0.3 || sum(voting) < 20
            t = t * 2;
            if t < 0.001
                do = 1; 
            end
        else
            do = 0;
        end        
    end
    score = score * dis;
    % v2
    
%     weight = voting ./ (voting_all);
%     weight(isnan(weight)) = 0;
%     weight(weight > 1) = 1;
%     score = weight * validMask;
end


% reprojection error.
function [inliers, score, H] = homogError(H, x, t, index, validMask)
    
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
    
    x_index = floor((index(inliers) - 1) / 10) + 1;
    y_index = index(inliers) - x_index * 10;
    dis = (max(x_index) - min(x_index) + 1) * (max(y_index) - min(y_index) + 1) / 100;
    
    voting = histcounts([index(inliers); 1; 100], 100);  
    voting(1) = voting(1) - 1;
    voting(100) = voting(100) - 1;
    
    voting_all = histcounts([index; 1; 100], 100);  
    voting_all(1) = voting_all(1) - 1;
    voting_all(100) = voting_all(100) - 1;
    
    voting(voting < voting_all * 0.5) = 0;
    voting(voting > 0) = 1;
    score = sum(voting & validMask') / sum(voting | validMask');
    
%     weight = voting ./ (voting_all);
%     weight(isnan(weight)) = 0;
%     weight(weight > 1) = 1;
%     score = weight * validMask;
    score = score * dis;
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
