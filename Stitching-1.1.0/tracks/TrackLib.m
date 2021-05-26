classdef TrackLib < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        videoWidth;
        videoHeight;
        nTrack;
        list;
        wSize;
        head;
        tail;
        id;
        len;
        points;
        labels;
        nFrame;
        live;
        nLabel;
        nWindow;
        labelMask;
    end
    
    methods
        function obj = TrackLib(nFrame)
            obj.points = [];
            obj.nFrame = nFrame;
            obj.nTrack = 0;
            obj.live = [];
            obj.head = [];
            obj.len = [];
            obj.tail = [];
        end        
        
        function setWindowSize(obj, size)
            obj.wSize = size; 
        end
        
        function addPoints(obj, newPoints, frameIndex)
            nNew = size(newPoints, 1);
            obj.points = [obj.points; zeros(nNew, obj.nFrame, 2)];            
            obj.points(obj.nTrack + 1:obj.nTrack + nNew, frameIndex, :) = reshape(newPoints, [nNew, 1, 2]);            
            obj.live = [obj.live obj.nTrack + 1: obj.nTrack + nNew];            
            obj.head = [obj.head; ones(nNew, 1) * frameIndex];
            obj.tail = [obj.tail; ones(nNew, 1) * frameIndex];
            obj.nTrack = obj.nTrack + nNew;
            obj.len = [obj.len; ones(nNew, 1)];
        end
        
        function updatePoints(obj, uPoints, frameIndex)
            obj.points(obj.live, frameIndex, :)  = reshape(uPoints, [length(obj.live), 1, 2]);
            obj.len(obj.live) = obj.len(obj.live) + 1;
        end
        
        function endPoints(obj, validity, frameIndex)
            obj.tail(obj.live(validity == false)) = frameIndex - 1;            
            obj.live(validity == false) = [];
        end
        
        function [M2, id] = getM(obj, index)
            left = (index - 1) * obj.wSize / 2 + 1;
            right = min(left + obj.wSize - 1, obj.nFrame);
            mid = (left + right) / 2.0;
            sub = obj.points(:, left:right, :) ;
            valid = 1:obj.nTrack;
            valid = valid(obj.head < mid & obj.tail > mid & obj.len > obj.wSize * 0.5);
            sub = sub(valid, :, :);
            id = valid;
            M2 = permute(sub, [3 1 2]);
        end
        
        function [f1, f2] = getF(obj, frameIndex)
%             selected = obj.head <= frameIndex & obj.tail > frameIndex & (obj.labels(:, frameIndex, 1) ==1  | obj.labels(:, frameIndex, 2) == 1);
            selected = obj.head <= frameIndex & obj.tail > frameIndex & obj.len > 4;
            f1 = squeeze(obj.points(selected, frameIndex, :));
            f2 = squeeze(obj.points(selected, frameIndex + 1, :));
            [~,inliners] = EstimateHomographyByRANSAC(f1',f2', 0.005);
            f1 = f1(inliners,:);
            f2 = f2(inliners,:);
        end
        function [f1, f2] = getAllF(obj, frameIndex)
%             selected = obj.head <= frameIndex & obj.tail > frameIndex & (obj.labels(:, frameIndex, 1) ==1  | obj.labels(:, frameIndex, 2) == 1);
            selected = obj.head <= frameIndex & obj.tail > frameIndex;
            f1 = squeeze(obj.points(selected, frameIndex, :));
            f2 = squeeze(obj.points(selected, frameIndex + 1, :));           
        end

        function [f1, f2] = getGoodF(obj, frameIndex, backList)
             selected = obj.head <= frameIndex & obj.tail > frameIndex & backList == 1;
             f1 = squeeze(obj.points(selected, frameIndex, :));
             f2 = squeeze(obj.points(selected, frameIndex + 1, :));
             
             [~,inliners] = EstimateHomographyByRANSAC(f1',f2', 0.001);
            f1 = f1(inliners,:);
            f2 = f2(inliners,:);
        end
        
        function addLabel(obj, wSize)
            epsilon = 0.1:0.1:0.4;
            obj.setWindowSize(wSize);
            obj.nWindow = ceil(obj.nFrame / (wSize / 2.0)) - 1;
            obj.labels = zeros(size(obj.points));
            nLb = 1;
            for windowIndex = 1:obj.nWindow
                tic;
                fprintf('window: %d - %d', (windowIndex - 1) * wSize / 2 + 1, (windowIndex - 1) * wSize / 2 + wSize);
                left = (windowIndex - 1) * wSize / 2 + 1;
                [rawData, trackID] = obj.getM(windowIndex);
                mask = (rawData > 0);
                mask = double(squeeze(mask(1,:,:) & mask(2,:,:))); 
                fprintf('\t %d tracks\n', size(rawData, 2));
                processedData = process_sequence(rawData, 'sparse', 'incomplete', mask);
                result = try_sequence('oc1R2RC', processedData, epsilon);
                computedLabels = find_best_segmentation(result, processedData, 100, epsilon);
                tempNlabel = max(computedLabels);
                if tempNlabel > nLb
                    nLb = tempNlabel; 
                end
                for label = 1:tempNlabel
                    obj.labels(trackID(computedLabels == label), left:left + wSize / 2 - 1, 2) = label;
                    obj.labels(trackID(computedLabels == label), left + wSize / 2:left + wSize - 1, 1) = label;
                end
                toc;
            end
            obj.nLabel = nLb;
        end
        
        function addFakeLabel(obj, wSize)
            obj.setWindowSize(wSize);
            obj.nWindow = ceil(obj.nFrame / (wSize / 2.0)) - 1;
            obj.labels = zeros(size(obj.points));
            obj.nLabel = 1;
            for windowIndex = 1:obj.nWindow
                tic;
                fprintf('window: %d - %d', (windowIndex - 1) * wSize / 2 + 1, (windowIndex - 1) * wSize / 2 + wSize);
                left = (windowIndex - 1) * wSize / 2 + 1;
                [~, trackID] = obj.getM(windowIndex);
                obj.labels(trackID, left:left + wSize / 2 - 1, 2) = 1;
                obj.labels(trackID, left + wSize / 2:left + wSize - 1, 1) = 1;
                toc;
            end
        end
        
        function refineLabel(obj, merge)
            nWindow = ceil(obj.nFrame / (obj.wSize / 2.0)) - 1; 
            nRound = 1;
            for roundIndex = 1:nRound
                for wIndex = 1:nWindow
                     
                end
            end
        end
        
        function LabelMask(obj)
            obj.labelMask = zeros(obj.nFrame, 2, obj.nLabel, obj.videoHeight, obj.videoWidth);
            for frameIndex = 1:obj.nFrame                
                hasPoint = obj.points(:, frameIndex, 1) ~= 0 | obj.points(:, frameIndex, 2) ~= 0;
                pointList = squeeze(tracks.points(hasPoint, frameIndex, :));
                labelList = squeeze(tracks.labels(hasPoint, frameIndex, :)); 
                for labelIndex = 1:obj.nLabel
                    mask = zeros(obj.videoHeight, obj.videoWidth);
                    mask((pointList(labelList(:, 1) == labelIndex, 2) - 1) * obj.videoWidth + pointList(labelList(:, 1) == labelIndex, 1)) = 10;
                    mask = imgaussfilt(mask);
                    mask(mask>1) = 1;                    
                    obj.labelMask(frameIndex, 1, labelIndex, :, :) = mask;
                    imshow(mask);
                    mask = zeros(obj.videoHeight, obj.videoWidth);
                    mask((pointList(labelList(:, 2) == labelIndex, 2) - 1) * obj.videoWidth + pointList(labelList(:, 2) == labelIndex, 1)) = 10;
                    mask = imgaussfilt(mask);
                    mask(mask>1) = 1;
                    obj.labelMask(frameIndex, 2, labelIndex, :, :) = mask;
                    imshow(mask);
                end
            end
        end
        
        function mask = getLabelMask(obj, frameIndex, leftright)
            hasPoint = obj.points(:, frameIndex, 1) ~= 0 | obj.points(:, frameIndex, 2) ~= 0;
            pointList = squeeze(obj.points(hasPoint, frameIndex, :));
            labelList = squeeze(obj.labels(hasPoint, frameIndex, leftright));
            mask = zeros(obj.videoHeight, obj.videoWidth, obj.nLabel);
            t_nLabel = max(labelList);
            for labelIndex = 1:t_nLabel
                mask1 = zeros(obj.videoHeight / 2, obj.videoWidth / 2);
                mask1((round(pointList(labelList == labelIndex, 1) / 2) - 1) * obj.videoHeight / 2 + round(pointList(labelList == labelIndex, 2) / 2)) = 1;
                mask_G = gpuArray(mask1);
                mask_G = imgaussfilt(mask_G, 50);                          
                mask_G = imgaussfilt(mask_G, 50);
                mask(:, :, labelIndex) = gather(imresize(mask_G, 2));
            end
            mask = mask * 1000;
            maskOld = mask;
            for labelIndex = 1:t_nLabel
                mask(:, :, labelIndex, 1) = maskOld(:, :, labelIndex, 1) ./ (sum(maskOld(:, :, :, 1), 3)+ 1e-100);                
            end
            drawMap = 0;
            if drawMap == 1
                colors = [0 0 1; 0 1 0; 1 0 0; 0 1 1; 1 0 1; 1 1 0; 1 1 1];
                merge1 = zeros(obj.videoHeight, obj.videoWidth, 3);

%                 [~, maxMask1] = max(mask(:, :, :, 1), [], 3);
%                 [~, maxMask2] = max(mask(:, :, :, 2), [], 3);
% 
%                 mask(:, :, 1, 1) = mask(:, :, 1, 1) .* double(maxMask1 == 1);
%                 mask(:, :, 2, 1) = mask(:, :, 2, 1) .* double(maxMask1 == 2);
%                 mask(:, :, 3, 1) = mask(:, :, 3, 1) .* double(maxMask1 == 3);
%                 mask(:, :, 4, 1) = mask(:, :, 4, 1) .* double(maxMask1 == 4);
%                 mask(:, :, 5, 1) = mask(:, :, 5, 1) .* double(maxMask1 == 5);            


                for labelIndex = 1:t_nLabel
                    binaryMask = mask(:, :, labelIndex);
                    submerge1 = repmat(binaryMask, [1 1 3]);
                    for chn = 1:3
                        submerge1(:, :, chn) = submerge1(:, :, chn) * colors(labelIndex, chn);
                    end
                    
%                     figure(1);imshow(submerge1);
%                     figure(2);imshow(submerge2);
                    merge1 = merge1 + submerge1;            
                end

                figure(3);imshow(merge1);
            end
        end
        function mask = getBackMask(obj, frameIndex, backList)
            backPoints = obj.head <= frameIndex & obj.tail >= frameIndex & backList == 1;
            forePoints = obj.head <= frameIndex & obj.tail >= frameIndex & backList == -1;
            mask1 = zeros(obj.videoHeight, obj.videoWidth);
            mask1((round(obj.points(backPoints, frameIndex, 1)) - 1) * obj.videoHeight + round(obj.points(backPoints, frameIndex, 2))) = 1;
            mask2 = zeros(obj.videoHeight, obj.videoWidth);
            mask2((round(obj.points(forePoints, frameIndex, 1)) - 1) * obj.videoHeight + round(obj.points(forePoints, frameIndex, 2))) = 1;
            mask1 = imgaussfilt(mask1, 50);
            mask1 = mask1 + imgaussfilt(mask1, 50);
            mask2 = imgaussfilt(mask2, 50);
            mask2 = mask2 + imgaussfilt(mask2, 50);
            mask = (mask1) ./ (mask1 + mask2);            
        end
        function mask = getForeMask(obj, frameIndex, backList)
            backPoints = obj.head <= frameIndex & obj.tail >= frameIndex & backList == 1;
            forePoints = obj.head <= frameIndex & obj.tail >= frameIndex & backList == -1;
            mask1 = zeros(obj.videoHeight, obj.videoWidth);
            mask1((round(obj.points(backPoints, frameIndex, 1)) - 1) * obj.videoHeight + round(obj.points(backPoints, frameIndex, 2))) = 1;
            mask2 = zeros(obj.videoHeight, obj.videoWidth);
            mask2((round(obj.points(forePoints, frameIndex, 1)) - 1) * obj.videoHeight + round(obj.points(forePoints, frameIndex, 2))) = 1;
            mask1 = imgaussfilt(mask1, 50);            
            mask2 = imgaussfilt(mask2, 50);            
            mask = (mask2) ./ (mask1 + mask2);            
            mask(isnan(mask)) = 0;
        end
        
        
    end
    
end

