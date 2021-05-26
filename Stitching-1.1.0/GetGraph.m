function [Path, graph, goodA, goodB] = GetGraph( trackA, trackB, CP, ppf, alpha, beta)
    nWindow = trackA.nWindow;
    nLabel = trackA.nLabel;
    videoH = trackA.videoHeight;
    videoW = trackA.videoWidth;
    graph = Graph(nWindow, nLabel);
    wSize = trackA.wSize;
    reGraph = 1;
    if reGraph == 0
        load('graph.mat') ;
    else        
        for windowIndex = 1:nWindow
            fprintf('\nwindow: %d\n', windowIndex);
            left = (windowIndex - 1) * wSize / 2 + 1;
            right= (windowIndex + 1) * wSize / 2; 
            nU_current = max(max(max(trackA.labels(:, left:left + wSize / 2 - 1, 2))), max(max(trackA.labels(:, left + wSize / 2:right, 1))));
            
            nV_current = max(max(max(trackB.labels(:, left:left + wSize / 2 - 1, 2))), max(max(trackB.labels(:, left + wSize / 2:right, 1))));
            
            
            if windowIndex < nWindow
                nU_next = max(max(max(trackA.labels(:, left + wSize / 2:left + wSize - 1, 2))), max(max(trackA.labels(:, left + wSize:right + wSize / 2, 1))));
                nV_next = max(max(max(trackB.labels(:, left + wSize / 2:left + wSize - 1, 2))), max(max(trackB.labels(:, left + wSize:right + wSize / 2, 1))));
    %             fprintf('inter-window\n');
                % inter-window
                left = windowIndex * wSize / 2 + 1;
                right= (windowIndex + 1) * wSize / 2;                
                for video = 1:2
                    if video == 1
                        track = trackA;
                        t_nLabel1 = nU_current;
                        t_nLabel2 = nU_next;
                    else
                        track = trackB;
                        t_nLabel1 = nV_current;
                        t_nLabel2 = nV_next;
                    end
                    for u = 1:t_nLabel1
                        for v = 1:t_nLabel2                        
                            % u in window windowIndex, v in window (windowIndex + 1) 
                            common = track.labels(:, left, 1) == u & track.labels(:, left, 2) == v;
                            common_ = track.labels(:, right, 1) == u & track.labels(:, right, 2) == v;
                            if common ~= common_
                                error('?'); % these two should be the same
                            end 
                            nCommon = sum(common);

                            if nCommon < 10
                                edgeWeight = 1e8;
                                graph.setEdgeWeight(nodeIndex(video, windowIndex, u, nLabel, nWindow), ...
                                    nodeIndex(video, windowIndex + 1, v, nLabel, nWindow), edgeWeight);
                                continue;
                            end

                            % compute the trajectory matrix
                            W2 = track.points(common, left:right, :);                    
                            a1 = W2(:, :, 1); a2 = W2(:, :, 2);
                            a1 = a1 ./ videoW; a2 = a2 ./ videoH;
                            a1a2 = [reshape(a1', [nCommon * wSize / 2 1]) reshape(a2', [nCommon * wSize / 2 1])];
                            W = reshape(a1a2', [wSize nCommon]);

                            fprintf('%d -> %d nCommon = %d, Rank = %d\n', u, v, nCommon, rank(W, 0.05));

                            % compute the edge weight
                            edgeWeight = exp( - alpha * nCommon) * rank(W, 0.05) / 4;
                            graph.setEdgeWeight(nodeIndex(video, windowIndex, u, nLabel, nWindow), ...
                                nodeIndex(video, windowIndex + 1, v, nLabel, nWindow), edgeWeight);                        
                        end
                    end
                end
            end
            % intra-window
    %         fprintf('intra-window\n=>');
            % 1. scan all the frames inside the window. 
            left = (windowIndex - 1) * wSize / 2 + 1;
            right= (windowIndex + 1) * wSize / 2; 
            edges = zeros(nLabel);
%             Mu = zeros(wSize, nLabel);
%             Mv = zeros(wSize, nLabel);
            deltaEdges = zeros(wSize, nLabel, nLabel);
            for fIndex = 1:wSize
                frameIndex = left + fIndex - 1;
                if (frameIndex > trackA.nFrame || frameIndex > trackB.nFrame) 
                   break;
                end
                labelMasks = zeros(videoH, videoW, nLabel, 2);
                fprintf('%4d', frameIndex);
                if frameIndex < left + wSize / 2
                    labelMasks(:, :, :, 1) = trackA.getLabelMask(frameIndex, 2);
                    labelMasks(:, :, :, 2) = trackB.getLabelMask(frameIndex, 2);
                else
                    labelMasks(:, :, :, 1) = trackA.getLabelMask(frameIndex, 1);
                    labelMasks(:, :, :, 2) = trackB.getLabelMask(frameIndex, 1);
                end

                pa = squeeze(round(CP(frameIndex, 1:ppf(frameIndex), 1:2))); % N x 2
                pb = squeeze(round(CP(frameIndex, 1:ppf(frameIndex), 3:4)));
                valid = pa(:, 1)>0 & pa(:, 1)<=videoW & pb(:, 1)> 0 & pb(:, 1)<=videoW ...
                    & pa(:, 2)>0 & pa(:, 2)<=videoH & pb(:, 2)> 0 & pb(:, 2)<=videoH;
                pa = pa(valid, :);
                pb = pb(valid, :);                

                for u = 1:nU_current
                    for v = 1:nV_current
                        % u in video A, v in video B
                        % 2. compute inter-video edge weight
                        uMask = squeeze(labelMasks(:, :, u, 1));
                        vMask = squeeze(labelMasks(:, :, v, 2));
                        Mu = uMask(videoH * (pa(:, 1) - 1) + pa(:, 2));
                        Mv = vMask(videoH * (pb(:, 1) - 1) + pb(:, 2));
                        deltaEdges(fIndex, u, v) = Mu' * Mv / wSize;
%                         edges(u, v) = edges(u, v) + Mu' * Mv / wSize;
                    end
                end
            end            
            for frameIndex = left:right
                for u = 1:nU_current
                    for v = 1:nV_current
                        edges(u, v) = edges(u, v) + deltaEdges(frameIndex - left + 1, u, v);
                    end                    
                end
            end

            for u = 1:nLabel
                for v = 1:nLabel
                    edgeWeight = exp( - beta * edges(u, v));
                    graph.setEdgeWeight(nodeIndex(1, windowIndex, u, nLabel, nWindow), ...
                        nodeIndex(2, windowIndex, v, nLabel, nWindow), edgeWeight);
                end
            end   
        end
        save('graph.mat', 'graph');
    end
%     graph.computeCoupledShortestPathfrom0tok();
    graph.computeAll();

    Path = graph.getShortestPath();
%     graph.computeCoupledShortestPathfromkto0();
    goodA = zeros(nWindow, nLabel);
    goodB = zeros(nWindow, nLabel);
    Path.dis
    for windowIndex = 1:nWindow
        fprintf('Window %2d\n', windowIndex);
        goodA(windowIndex, Path.video1Relative(windowIndex) + 1) = 1;
        fprintf('%d\t%f*\n', Path.video1Relative(windowIndex) + 1, Path.dis1);
       
        for testLabel = 1:nLabel                
            node = nodeIndex(1, windowIndex, testLabel, nLabel, nWindow);    
            Path1 = graph.getShortestPathAcrossSpecificNode(node);
            dis = Path1.dis;      
            fprintf('%d\t%f\t', testLabel, dis);
            if dis < max(2 * Path.dis1, 2)
                goodA(windowIndex, testLabel) = 1;                
                fprintf('[');
                [left, right] = getNearest(graph, 1, windowIndex, testLabel, nLabel, nWindow);
                if windowIndex > 1 && left > 0
                    goodA(windowIndex - 1, left) = 1;
                    fprintf('%2d', left);
                end
                if windowIndex < nWindow && right > 0
                    goodA(windowIndex + 1, right) = 1; 
                    fprintf('%2d', right);
                end 
                fprintf(']');
            end            
        end
        
        fprintf('\n');
        goodB(windowIndex, Path.video2Relative(windowIndex) + 1) = 1;
        fprintf('%d\t%f*\n', Path.video2Relative(windowIndex) + 1, Path.dis2);
        
        for testLabel = 1:nLabel
            node = nodeIndex(2, windowIndex, testLabel, nLabel, nWindow);    
            Path2 = graph.getShortestPathAcrossSpecificNode(node);
            dis = Path2.dis;                
            fprintf('%d\t%f\t', testLabel, dis);
            if dis < max(2 * Path.dis2, 2)
                goodB(windowIndex, testLabel) = 1;   
                fprintf('[');
                [left, right] = getNearest(graph, 2, windowIndex, testLabel, nLabel, nWindow);
                if windowIndex > 1 && left > 0;
                    fprintf('%2d', left);
                    goodB(windowIndex - 1, left) = 1;
                end
                if windowIndex < nWindow && right > 0;
                    fprintf('%2d', right);
                    goodB(windowIndex + 1, right) = 1; 
                end
                fprintf(']');
            end
        end
        fprintf('\n');
        
        
%         graph_new = copy(graph);        
%         graph_new.ignoreNode(nodeIndex(1, windowIndex, Path.video1Relative(windowIndex) + 1, nLabel, nWindow));
%         graph_new.computeCoupledShortestPathfrom0tok();
%         path_new = graph_new.getShortestPath();
%         wIndex = 1:nWindow;
%         path_new.dis
%         
%         
%         
%         goodB(windowIndex, Path.video2Relative(windowIndex) + 1) = 1;
%         graph_new = copy(graph);
%         graph_new.ignoreNode(nodeIndex(2, windowIndex, Path.video2Relative(windowIndex) + 1, nLabel, nWindow));
%         graph_new.computeCoupledShortestPathfrom0tok();
%         path_new = graph_new.getShortestPath();
%         path_new.dis
%         wIndex = 1:nWindow;
%         if path_new.dis < 1e30 * Path.dis
%             goodA(nWindow * path_new.video1Relative(wIndex) + wIndex) = 1; 
%             goodB(nWindow * path_new.video2Relative(wIndex) + wIndex) = 1;
%         end
    end
    
end

function index = nodeIndex(v, w, u, nL, nW)
    index = (v-1)*nL*nW + (w-1)*nL + u;
end



function new = copy(this)
    % Instantiate new object of the same class.
    new = feval(class(this));

    % Copy all non-hidden properties.
    p = properties(this);
    for i = 1:length(p)
        new.(p{i}) = this.(p{i});
    end
end

function dis = getShortestDisBy(graph, videoIndex, windowIndex, labelIndex, nLabel, nWindow)
    distances = zeros(nLabel, 1);
    node = nodeIndex(videoIndex, windowIndex, labelIndex, nLabel, nWindow);    
    for labelIndex = 1:nLabel
        nodeTest = nodeIndex(3 - videoIndex, windowIndex, labelIndex, nLabel, nWindow);   
        if videoIndex == 1
            path_new = graph.getShortestPathAcrossSpecificNodes(node, nodeTest);
        else
            path_new = graph.getShortestPathAcrossSpecificNodes(nodeTest, node);
        end
        distances(labelIndex) = path_new.dis;
    end
    dis = min(distances);
end

function [left, right] = getNearest(graph, videoIndex, windowIndex, labelIndex, nLabel, nWindow)    
    node = nodeIndex(videoIndex, windowIndex, labelIndex, nLabel, nWindow);
    left = 0;
    right = 0;
    if windowIndex > 1
        [edge, nearest] = min(squeeze(graph.adjacentMatrix(node, node - labelIndex - nLabel + 1:node - labelIndex)));
        if edge < 1e8
            left = nearest;
        end
    end
    if windowIndex < nWindow
        [edge, nearest] = min(squeeze(graph.adjacentMatrix(node, node - labelIndex + nLabel + 1:node - labelIndex + nLabel + nLabel)));
        if edge < 1e8
            right = nearest;
        end
    end
end
