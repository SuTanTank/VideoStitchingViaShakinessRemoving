
function [M, inliers] = RANSAC2(x, fittingfn, distfn, degenfn, s, t, index, validMask)

    maxTrials = 1000;
    maxDataTrials = 100;
    
    [rows, npts] = size(x);
    
    % Desired probability of choosing at least one sample free from outliers (probably should be a parameter)
    p = 0.99; 

    bestM = NaN;      
    trialcount = 0;
    bestscore =  0;
    N = 1;            
    
    while N > trialcount
        
        % Select at random s datapoints to form a trial model, M.
        degenerate = 1;
        count = 1;
        while degenerate
            % Generate s random indicies in the range 1..npts
            ind = zeros(s, 1);
%             ind = randsample(npts, s);
            for i = 1:s
                pick = randsample(npts, 1);
                while validMask(index(pick)) == 0 || any(ind == pick) || any(index(ind(ind > 0)) == index(pick))
                    pick = randsample(npts, 1);
                end
                ind(i) = pick;
            end

            % Test that these points are not a degenerate configuration.
            degenerate = feval(degenfn, x(:,ind));
            
            if ~degenerate
                M = feval(fittingfn, x(:,ind));
                if isempty(M)
                    degenerate = 1;
                end
            end
            
            % Safeguard against being stuck in this loop forever
            count = count + 1;
            if count > maxDataTrials
                disp('Unable to select a nondegenerate data set');
                break
            end
        end
        
%         if scoreX1 ~= 0
%             l = min([scoreX1(ind) scoreX2(ind)]);
%             t = (l + 0.6) * t;
%         end
        
        [inliers, score, M] = feval(distfn, M, x, t, index, validMask);
        
        % Find the number of inliers to this model.
        ninliers = length(inliers);
        
        if score >= bestscore   
            bestscore = score; 
            bestinliers = inliers;
            bestM = M;
            
            % Update estimate of N
            fracinliers =  ninliers/npts;
            pNoOutliers = 1 -  fracinliers^s;
            pNoOutliers = max(eps, pNoOutliers);  
            pNoOutliers = min(1-eps, pNoOutliers);
            N = log(1-p)/log(pNoOutliers);
        end
        
        trialcount = trialcount+1;
        
        % Safeguard against being stuck in this loop forever
        if trialcount > maxTrials
            fprintf('ransac reached the maximum number of %d trials\n', maxTrials);
            break
        end
        
        %fprintf('Iter: %d \n', trialcount);
    end
    
    if ~isnan(bestM)
        M = bestM;
        inliers = bestinliers;
    else
        M = [];
        inliers = [];
        disp('ransac was unable to find a useful solution');
    end
end
