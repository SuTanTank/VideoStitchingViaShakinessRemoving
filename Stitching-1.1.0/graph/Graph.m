classdef Graph < handle
    %GRAPH Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        %% common variables
        adjacentMatrix;     % The adjacentMatrix which stores the distance between nodes
        clusterNum;         % The maximum number of clusters in one temporal window in one video
        windowNum;          % The number of the temporal windows        
        nodeNum;            % Intermediate value: The total number of nodes in the graph        
        
        %% common results
        A;                  % assist cells for forward shortest path: Node1Num * Node2Num
        B;                  % assist cells for backward shortest path: Node1Num * Node2Num        
        ShortestPath;       % shortest path
        ShortestPathReverse; % reverse shortest path
        
        %% individual results
        % assist variables for video1
        A1;
        B1;
        ShortestPath1;
        ShortestPathReverse1;
        
        % assist variables for video2
        A2;
        B2;
        ShortestPath2;
        ShortestPathReverse2;        
        
        %% common constant values
        maxDistance;        % Constant value: 1e10
        dftDistance;        % Constant value: 1e8
    end
    
    methods
         function obj = Graph(windowNum, clusterNum)
             if nargin > 2
                 error('Graph construcor too many inputs, require at most 2 optional inputs');
             end
             
             switch nargin
                 case 0
                     obj.windowNum = 24;
                     obj.clusterNum = 5;
                 case 1
                     obj.windowNum = windowNum;
                     obj.clusterNum = 5;
                 case 2
                     obj.windowNum = windowNum;
                     obj.clusterNum = clusterNum;
             end
             
             obj.maxDistance = 1e10; 
             obj.dftDistance = 1e6;    % 1e8
                          
             obj.nodeNum = 2 * obj.clusterNum * obj.windowNum;
             
             obj.A = cell(obj.nodeNum, obj.nodeNum);
             obj.B = cell(obj.nodeNum, obj.nodeNum);
             
             % We purposely set the default distance between two nodes to
             % be very large
             obj.adjacentMatrix = ones(obj.nodeNum, obj.nodeNum) * obj.dftDistance;
             
             obj.A1 = cell(obj.clusterNum * obj.windowNum, 1);   
             obj.B1 = cell(obj.clusterNum * obj.windowNum, 1);
             obj.A2 = cell(obj.clusterNum * obj.windowNum, 1);
             obj.B2 = cell(obj.clusterNum * obj.windowNum, 1);
         end
         
         function idx = getNodeIdx(obj, nodeVideoIdx, nodeWinIdx, nodeIdx)
             % GETNODEIDX Get the index of the node in the adjecent matrix
             % nodeVideoIdx: 0 or 1, indicating which video
             % nodeWinIdx: 0,1,...,obj.windowNum-1, indicating which
             % temporal window
             % nodeIdx: 0,1,...,obj.clusterNum -1, indicating which node
             % 
             % the order of the node is as follows:
             %      A      |      B
             % w1       w2 | w1       w2
             % node1w1A, node2w1A, node3w1A, node1w2A, node2w2A, node3w2A | 
             % node1w1B, node2w1B, node3w1B, node1w2B, node2w2B, node3w2B 
             
             if nodeVideoIdx > 1
                error('GetNodeIdx: nodeVideoIdx should be 0 or 1');
             end
             
             if nodeWinIdx > obj.windowNum - 1
                error('GetNodeIdx: nodeWinIdx should be 0,1,...,obj.windowNum-1');
             end
             
             if nodeIdx > obj.clusterNum - 1
                error('GetNodeIdx: nodeIdx should be 0,1,...,obj.clusterNum-1');
             end
             
             videoIdxOffset = nodeVideoIdx * obj.windowNum * obj.clusterNum;
             winIdxOffset = nodeWinIdx * obj.clusterNum;
             
             idx = videoIdxOffset + winIdxOffset + nodeIdx + 1;             
         end
         
         function ignoreRelativeNode(obj, nodeVideoIdx, nodeWinIdx, nodeIdx)
             % ignoreNode set all the edges connected to this node to be
             % very large
             
             idx = getNodeIdx(obj, nodeVideoIdx, nodeWinIdx, nodeIdx);
             
             obj.adjacentMatrix(idx, :) = obj.maxDistance;
             obj.adjacentMatrix(:, idx) = obj.maxDistance;
         end
         
         function ignoreNode(obj, idx)
             % ignoreNode set all the edges connected to this node to be
             % very large
             
             obj.adjacentMatrix(idx, :) = obj.maxDistance;
             obj.adjacentMatrix(:, idx) = obj.maxDistance;
         end

         function [videoIdx, winIdx, nodeIdx] = getRelativeIdx(obj, idx)
            %  getRelativeIdx This function converts absolute index to
            %  relative index
            idx = idx - 1;
            
            videoIdx = floor(idx/( obj.windowNum * obj.clusterNum));
            winIdx = floor((idx - videoIdx * obj.windowNum * obj.clusterNum)/ obj.clusterNum);
            nodeIdx = idx - videoIdx * obj.windowNum * obj.clusterNum - winIdx * obj.clusterNum;
         end
         
         function setEdgeWeightRelative(obj, node1VideoIdx, node1WinIdx, node1Idx, node2VideoIdx, node2WinIdx, node2Idx, weight)
            % setEdgeWeightRelative set the weights of edges between nodes, using
            % relative index, starts from 0
            idx1 = getNodeIdx(obj, node1VideoIdx, node1WinIdx, node1Idx);
            idx2 = getNodeIdx(obj, node2VideoIdx, node2WinIdx, node2Idx);
            
            obj.adjacentMatrix(idx1, idx2) = weight;
            obj.adjacentMatrix(idx2, idx1) = weight;
         end
         
         function weight = getEdgeWeightRelative(obj, node1VideoIdx, node1WinIdx, node1Idx, node2VideoIdx, node2WinIdx, node2Idx)
            % getEdgeWeightRelative get the weights of edges between nodes, using
            % relative index
            idx1 = getNodeIdx(obj, node1VideoIdx, node1WinIdx, node1Idx);
            idx2 = getNodeIdx(obj, node2VideoIdx, node2WinIdx, node2Idx);
            
            weight = obj.adjacentMatrix(idx1, idx2);
         end
         
         
         function setEdgeWeight(obj, node1Idx, node2Idx, weight)
             % setEdgeWeight set the weights of edges between nodes, using absolute index
             % starts from 1
            obj.adjacentMatrix(node1Idx, node2Idx) = weight;
            obj.adjacentMatrix(node2Idx, node1Idx) = weight;
         end
         
         function weight = getEdgeWeight(obj, node1Idx, node2Idx)
             % getEdgeWeight get the weights of edges between nodes, using absolute index
             weight = obj.adjacentMatrix(node1Idx, node2Idx);
         end
         
         function generateRNDWeights(obj)
             % GENERATERNNDWEIGHTS this function generates random weights
             % for each edge
             
             rng(0,'twister');
             
             % for the 0th video
             nodeVideoIdx = 0;
             for nodeWinIdx = 0:obj.windowNum - 2
                for nodeIdx = 0:obj.clusterNum - 1
                    idx1 = getNodeIdx(obj, nodeVideoIdx, nodeWinIdx, nodeIdx);
                    for nodeIdxNextWin = 0:obj.clusterNum -1
                        idx2 = getNodeIdx(obj, nodeVideoIdx, nodeWinIdx + 1, nodeIdxNextWin);
                        
                        weight = randi([1 10]);
                        
                        obj.adjacentMatrix(idx1, idx2) = weight;
                        obj.adjacentMatrix(idx2, idx1) = weight;
                    end
                end
             end
             
             % for the 1th video
             nodeVideoIdx = 1;
             for nodeWinIdx = 0:obj.windowNum - 2
                for nodeIdx = 0:obj.clusterNum - 1
                    idx1 = getNodeIdx(obj, nodeVideoIdx, nodeWinIdx, nodeIdx);
                    for nodeIdxNextWin = 0:obj.clusterNum -1
                        idx2 = getNodeIdx(obj, nodeVideoIdx, nodeWinIdx + 1, nodeIdxNextWin);
                        
                        weight = randi([1 10]);
                        
                        obj.adjacentMatrix(idx1, idx2) = weight;
                        obj.adjacentMatrix(idx2, idx1) = weight;
                    end
                end
             end
             
             % for the connections between two videos
             for nodeWinIdx = 0:obj.windowNum -1
                for node1Idx = 0:obj.clusterNum - 1
                    idx1 = getNodeIdx(obj, 0, nodeWinIdx, node1Idx);
                    for node2Idx = 0:obj.clusterNum - 1
                        idx2 = getNodeIdx(obj, 1, nodeWinIdx, node2Idx);
                        
                        weight = randi([1 10]);
                        
                        obj.adjacentMatrix(idx1, idx2) = weight;
                        obj.adjacentMatrix(idx2, idx1) = weight;
                    end
                end
             end 
         end
         
         
         function computeCoupledShortestPathfrom0tok(obj)
            % computeCoupledShortestPathfrom0tok we compute the shortest
            % path using DP in this function
            % dis(k,s,t)=min{dis(k-1, m, n) + w(k-1,m,k,s) + w(k-1,n,k,t) + w(k,s,k,t)}
            
            k = 0; 
            node1VideoIdx = 0;
            node2VideoIdx = 1;
            
            for s = 0:(obj.clusterNum - 1)
                idx1 = getNodeIdx(obj, node1VideoIdx, k, s);
                for t = 0:(obj.clusterNum - 1)
                    idx2 = getNodeIdx(obj, node2VideoIdx, k, t);
                    obj.A{idx1, idx2}.dis = obj.adjacentMatrix(idx1, idx2);
                    obj.A{idx1, idx2}.m = idx1;
                    obj.A{idx1, idx2}.n = idx2;
                end
            end
            
            for k = 1:(obj.windowNum - 1)
                for s = 0:(obj.clusterNum - 1)
                    idx1 = getNodeIdx(obj, node1VideoIdx, k, s);
                    for t = 0:(obj.clusterNum - 1)
                       idx2 = getNodeIdx(obj, node2VideoIdx, k, t); 
                       
                       obj.A{idx1, idx2}.dis = obj.maxDistance;     %ideally, we would like to set it the positive infinity
                       obj.A{idx1, idx2}.m = -1;
                       obj.A{idx1, idx2}.n = -1;
                       
                       
                       for m = 0:(obj.clusterNum - 1)
                           idx3 = getNodeIdx(obj, node1VideoIdx, k - 1, m);
                           for n = 0:(obj.clusterNum - 1)
                               idx4 = getNodeIdx(obj, node2VideoIdx, k - 1, n);
                               
                               tmp = obj.A{idx3, idx4}.dis + obj.adjacentMatrix(idx3, idx1) + obj.adjacentMatrix(idx4, idx2) + obj.adjacentMatrix(idx1, idx2);
                               if tmp < obj.A{idx1, idx2}.dis
                                   obj.A{idx1, idx2}.dis = tmp;
                                   obj.A{idx1, idx2}.m = idx3;
                                   obj.A{idx1, idx2}.n = idx4;
                               end
                               
                           end
                       end
                       
                    end %t
                end %s
            end %k  
         end
         
         function dis = computePathDistance(obj, pathNodeIdx)
            dis = 0;
            numOfNodes = numel(pathNodeIdx);
            for i = 1:(numOfNodes - 1)
                dis = dis + obj.adjacentMatrix(pathNodeIdx(i), pathNodeIdx(i + 1));
            end
        end
         
         function Path = getShortestPath(obj)
            % return the shortest paths
            % Path.dis: the length of the shortest path
            % Path.video1 = []: a list of vertices index from the first
            % window to the last window for video 1
            % Path.video2 = []: a list of vertices index from the first
            % window to the last window for video 2
            
            tmpIdx1 = 0;
            tmpIdx2 = 0;
            tmpDis = obj.maxDistance;
            
            node1VideoIdx = 0;
            node2VideoIdx = 1;
            
            k = obj.windowNum -1;
            
            for i = 0:(obj.clusterNum - 1)
                idx1 = getNodeIdx(obj, node1VideoIdx, k, i);
                for j = 0:(obj.clusterNum - 1)
                    idx2 = getNodeIdx(obj, node2VideoIdx, k, j);
                    
                    if obj.A{idx1, idx2}.dis < tmpDis
                        tmpDis = obj.A{idx1, idx2}.dis;
                        tmpIdx1 = idx1;
                        tmpIdx2 = idx2;
                    end
                end
            end
            
            Path.dis = tmpDis;
            Path.video1 = tmpIdx1;
            Path.video2 = tmpIdx2;
            
            [~, ~, tmpNode1Idx] = getRelativeIdx(obj, tmpIdx1);
            [~, ~, tmpNode2Idx] = getRelativeIdx(obj, tmpIdx2);
            
            Path.video1Relative = tmpNode1Idx;
            Path.video2Relative = tmpNode2Idx;
            
            idx1 = tmpIdx1;
            idx2 = tmpIdx2;
            
            for k = obj.windowNum:-1:2
                tmpIdx1 = obj.A{idx1, idx2}.m;
                tmpIdx2 = obj.A{idx1, idx2}.n;
                
                Path.video1 = [tmpIdx1, Path.video1];
                Path.video2 = [tmpIdx2, Path.video2];
                
                [~, ~, tmpNode1Idx] = getRelativeIdx(obj, tmpIdx1);
                [~, ~, tmpNode2Idx] = getRelativeIdx(obj, tmpIdx2);
                
                Path.video1Relative = [tmpNode1Idx, Path.video1Relative];
                Path.video2Relative = [tmpNode2Idx, Path.video2Relative];
                
                idx1 = tmpIdx1;
                idx2 = tmpIdx2;
            end
            
            % compute the distance for each video
            Path.dis1 = computePathDistance(obj, Path.video1);
            Path.dis2 = computePathDistance(obj, Path.video2);
            
            obj.ShortestPath = Path;
         end
         
         function Path = computePathValueusingPaths(obj, video1, video2)
             % compute a path value according to the selected path
                      
             node1VideoIdx = 0;
             node2VideoIdx = 1;
             
             idx1 = getNodeIdx(obj, node1VideoIdx, 0, video1(1));
             idx2 = getNodeIdx(obj, node2VideoIdx, 0, video2(1));
             
             Path.dis = obj.adjacentMatrix(idx1, idx2);
             Path.video1 = idx1;
             Path.video2 = idx2;
             
             Path.video1Relative = video1;
             Path.video2Relative = video2;
             
             for i1 = 1:(obj.windowNum - 1)
                 preIdx1 = idx1;
                 idx1 = getNodeIdx(obj, node1VideoIdx, i1, video1(i1 + 1));

               
                 preIdx2 = idx2;
                 idx2 = getNodeIdx(obj, node2VideoIdx, i1, video2(i1 + 1));

                 Path.dis = Path.dis + obj.adjacentMatrix(idx1, idx2) + obj.adjacentMatrix(idx1, preIdx1) + obj.adjacentMatrix(idx2, preIdx2);
                 Path.video1 = [Path.video1 idx1];
                 Path.video2 = [Path.video2 idx2];

             end 
                 
         end
         
         function Path = getShortestPathBF(obj)
            % getShortestPathBF we use brute force to compute a shortest
             % path, this function is for validation. Be careful of using
             % this method, as the running time will increase
             % exponentially. O(clusterNum ^ windowNum)
             
             Path.dis = obj.maxDistance;
             Path.video1 = zeros(obj.windowNum);
             Path.video2 = zeros(obj.windowNum);
             
             numOfCases = int64(obj.clusterNum ^ (2*obj.windowNum));
             
             for caseIdx = 0:(numOfCases - 1)
                 caseCode = dec2base(caseIdx, obj.clusterNum, 2*obj.windowNum);
                 nodeIndices = caseCode - '0';
                 
                 tmpVideo1 = nodeIndices(1:obj.windowNum);
                 tmpVideo2 = nodeIndices(obj.windowNum + 1:2*obj.windowNum);
                                  
%                  disp(caseCode);

                 tmpPath = computePathValueusingPaths(obj, tmpVideo1, tmpVideo2);
                 if tmpPath.dis < Path.dis
                    Path = tmpPath;
                 end
                 
             end
             
         end
         
         
         function Path = getShortestPathBF2(obj)
             % getShortestPathBF we use brute force to compute a shortest
             % path, this function is for validation. Be careful of using
             % this method, as the running time will increase
             % exponentially. O(clusterNum ^ windowNum)
             % This function works only for windowNum = 2.
             
             Path.dis = obj.maxDistance;
             Path.video1 = zeros(obj.windowNum);
             Path.video2 = zeros(obj.windowNum);
             
             for i1 = 0:(obj.clusterNum - 1)
                 idx1 = getNodeIdx(obj, 0, 0, i1); 
                for i2 = 0:(obj.clusterNum - 1)
                    idx2 = getNodeIdx(obj, 0, 1, i2);
                    for i3 = 0:(obj.clusterNum - 1)
                        idx3 = getNodeIdx(obj, 1, 0, i3);
                        for i4 = 0:(obj.clusterNum - 1)
                            idx4 = getNodeIdx(obj, 1, 1, i4);
                            
                            tmpDis = obj.adjacentMatrix(idx1, idx2) + obj.adjacentMatrix(idx1, idx3) + obj.adjacentMatrix(idx3, idx4) + obj.adjacentMatrix(idx2, idx4);
                            
                            if tmpDis < Path.dis
                                Path.dis = tmpDis;
                                Path.video1 = [idx1, idx2];
                                Path.video2 = [idx3, idx4];
                            end
                        end
                    end
                end
             end
         end
         
        function computeCoupledShortestPathfromkto0(obj)
            % computeCoupledShortestPathfrom0tok we compute the shortest
            % path using DP in this function
            % dis(k,s,t)=min{dis(k-1, m, n) + w(k-1,m,k,s) + w(k-1,n,k,t) + w(k,s,k,t)}
            
            k = obj.windowNum; 
            node1VideoIdx = 0;
            node2VideoIdx = 1;
            
            for s = 0:(obj.clusterNum - 1)
                idx1 = getNodeIdx(obj, node1VideoIdx, k - 1, s);
                for t = 0:(obj.clusterNum - 1)
                    idx2 = getNodeIdx(obj, node2VideoIdx, k - 1, t);
                    obj.B{idx1, idx2}.dis = obj.adjacentMatrix(idx1, idx2);
                    obj.B{idx1, idx2}.m = idx1;
                    obj.B{idx1, idx2}.n = idx2;
                end
            end
            
            for k = (obj.windowNum - 1):-1:1
                for s = 0:(obj.clusterNum - 1)
                    idx1 = getNodeIdx(obj, node1VideoIdx, k - 1, s);
                    for t = 0:(obj.clusterNum - 1)
                       idx2 = getNodeIdx(obj, node2VideoIdx, k - 1, t); 
                       
                       obj.B{idx1, idx2}.dis = obj.maxDistance;     %ideally, we would like to set it the positive infinity
                       obj.B{idx1, idx2}.m = -1;
                       obj.B{idx1, idx2}.n = -1;
                                             
                       for m = 0:(obj.clusterNum - 1)
                           idx3 = getNodeIdx(obj, node1VideoIdx, k, m);
                           for n = 0:(obj.clusterNum - 1)
                               idx4 = getNodeIdx(obj, node2VideoIdx, k, n);
                               
                               tmp = obj.B{idx3, idx4}.dis + obj.adjacentMatrix(idx3, idx1) + obj.adjacentMatrix(idx4, idx2) + obj.adjacentMatrix(idx1, idx2);
                               if tmp < obj.B{idx1, idx2}.dis
                                   obj.B{idx1, idx2}.dis = tmp;
                                   obj.B{idx1, idx2}.m = idx3;
                                   obj.B{idx1, idx2}.n = idx4;
                               end
                               
                           end
                       end
                       
                    end %t
                end %s
            end %k  
        end
         
        function Path = getShortestPathfromkto0(obj)
            % return the shortest paths
            % Path.dis: the length of the shortest path
            % Path.video1 = []: a list of vertices index from the first
            % window to the last window for video 1
            % Path.video2 = []: a list of vertices index from the first
            % window to the last window for video 2
            
            tmpIdx1 = 0;
            tmpIdx2 = 0;
            tmpDis = obj.maxDistance;
            
            node1VideoIdx = 0;
            node2VideoIdx = 1;
            
            k = 1;
            
            for i = 0:(obj.clusterNum - 1)
                idx1 = getNodeIdx(obj, node1VideoIdx, k - 1, i);
                for j = 0:(obj.clusterNum - 1)
                    idx2 = getNodeIdx(obj, node2VideoIdx, k - 1, j);
                    
                    if obj.B{idx1, idx2}.dis < tmpDis
                        tmpDis = obj.B{idx1, idx2}.dis;
                        tmpIdx1 = idx1;
                        tmpIdx2 = idx2;
                    end
                end
            end
            
            Path.dis = tmpDis;
            Path.video1 = tmpIdx1;
            Path.video2 = tmpIdx2;
            
            [~, ~, tmpNode1Idx] = getRelativeIdx(obj, tmpIdx1);
            [~, ~, tmpNode2Idx] = getRelativeIdx(obj, tmpIdx2);
            
            Path.video1Relative = tmpNode1Idx;
            Path.video2Relative = tmpNode2Idx;
            
            idx1 = tmpIdx1;
            idx2 = tmpIdx2;
            
            for k = 1:obj.windowNum - 1
                tmpIdx1 = obj.B{idx1, idx2}.m;
                tmpIdx2 = obj.B{idx1, idx2}.n;
                
                Path.video1 = [tmpIdx1, Path.video1];
                Path.video2 = [tmpIdx2, Path.video2];
                
                [~, ~, tmpNode1Idx] = getRelativeIdx(obj, tmpIdx1);
                [~, ~, tmpNode2Idx] = getRelativeIdx(obj, tmpIdx2);
                
                Path.video1Relative = [tmpNode1Idx, Path.video1Relative];
                Path.video2Relative = [tmpNode2Idx, Path.video2Relative];
                
                idx1 = tmpIdx1;
                idx2 = tmpIdx2;
            end
            
            Path.video1 = flip(Path.video1);
            Path.video2 = flip(Path.video2);
            Path.video1Relative = flip(Path.video1Relative);
            Path.video2Relative = flip(Path.video2Relative);
            
            obj.ShortestPathReverse = Path;
        end
         
        function Path = getShortestPathAcrossSpecificNodes(obj, idx1, idx2)
            % getShortestPathAcrossSpecificNodes use absolute index, let
            % idx1 be a node from video1, and idx2 be a node from video2
            
            [tmpNode1VideoIdx, tmpNode1WindowIdx, tmpNode1Idx] = getRelativeIdx(obj, idx1);
            [tmpNode2VideoIdx, tmpNode2WindowIdx, tmpNode2Idx] = getRelativeIdx(obj, idx2);
            
            if tmpNode1WindowIdx ~= tmpNode2WindowIdx
                error('Window indices of the two nodes do no match!');
            end
            
            if tmpNode1VideoIdx == tmpNode2VideoIdx
                error('The nodes should come from different videos!');
            end
            
            Path.dis = obj.A{idx1, idx2}.dis + obj.B{idx1, idx2}.dis - obj.adjacentMatrix(idx1, idx2);
            Path.video1 = idx1;
            Path.video2 = idx2;
            
            %% construct a shortest path
            
            % part 1
            curIdx1 = idx1;
            curIdx2 = idx2;            
            
            Path.video1Relative = tmpNode1Idx;
            Path.video2Relative = tmpNode2Idx;
            
            for k = tmpNode1WindowIdx:-1:1
                tmpIdx1 = obj.A{curIdx1, curIdx2}.m;
                tmpIdx2 = obj.A{curIdx1, curIdx2}.n;
                
                Path.video1 = [tmpIdx1, Path.video1];
                Path.video2 = [tmpIdx2, Path.video2];
                
                [~, ~, tmpNode1Idx] = getRelativeIdx(obj, tmpIdx1);
                [~, ~, tmpNode2Idx] = getRelativeIdx(obj, tmpIdx2);
                
                Path.video1Relative = [tmpNode1Idx, Path.video1Relative];
                Path.video2Relative = [tmpNode2Idx, Path.video2Relative];
                
                curIdx1 = tmpIdx1;
                curIdx2 = tmpIdx2;
            end
            
            % part 2
            curIdx1 = idx1;
            curIdx2 = idx2; 
            for k = tmpNode1WindowIdx:obj.windowNum - 2
                tmpIdx1 = obj.B{curIdx1, curIdx2}.m;
                tmpIdx2 = obj.B{curIdx1, curIdx2}.n;
                
                Path.video1 = [Path.video1, tmpIdx1];
                Path.video2 = [Path.video2, tmpIdx2];
                
                [~, ~, tmpNode1Idx] = getRelativeIdx(obj, tmpIdx1);
                [~, ~, tmpNode2Idx] = getRelativeIdx(obj, tmpIdx2);
                
                Path.video1Relative = [Path.video1Relative, tmpNode1Idx];
                Path.video2Relative = [Path.video2Relative, tmpNode2Idx];
                
                curIdx1 = tmpIdx1;
                curIdx2 = tmpIdx2;
            end
        end
        
        function Path = getShortestPathAcrossSpecificNodesRelative(obj, node1VideoIdx, node1WinIdx, node1Idx, node2VideoIdx, node2WinIdx, node2Idx)
            % getShortestPathAcrossSpecificNodes use relative index, let
            % idx1 be a node from video1, and idx2 be a node from video2
            idx1 = getNodeIdx(obj, node1VideoIdx, node1WinIdx, node1Idx);
            idx2 = getNodeIdx(obj, node2VideoIdx, node2WinIdx, node2Idx);
            
            Path = getShortestPathAcrossSpecificNodes(obj, idx1, idx2);
        end
        
        %% Deal with single video        
        function computeShortestPathForVideo1From0tok(obj)
            % compute shortest path according to siggraph asia 2016 ZHANG
            % Fanglve
            
            k = 0;
            node1VideoIdx = 0;
            
            for s = 0:(obj.clusterNum - 1)
                idx1 = getNodeIdx(obj, node1VideoIdx, k, s);
                obj.A1{idx1}.dis = 0;
                obj.A1{idx1}.m = -1;
            end
            
            for k = 1:(obj.windowNum - 1)
                for s = 0:(obj.clusterNum - 1)
                    idx1 = getNodeIdx(obj, node1VideoIdx, k, s);
                    obj.A1{idx1}.dis = obj.maxDistance;
                    obj.A1{idx1}.m = -1;
                    
                    for m = 0:(obj.clusterNum - 1)
                        idx2 = getNodeIdx(obj, node1VideoIdx, k - 1, m);
                        tmp = obj.A1{idx2}.dis + obj.adjacentMatrix(idx1, idx2);
                        
                        if tmp < obj.A1{idx1}.dis
                            obj.A1{idx1}.dis = tmp;
                            obj.A1{idx1}.m = idx2;
                        end
                    end
                end
            end
            
        end
        
        function computeShortestPathForVideo1Fromkto0(obj)
            k = obj.windowNum;
            node1VideoIdx = 0;
            
            for s = 0:(obj.clusterNum - 1)
                idx1 = getNodeIdx(obj, node1VideoIdx, k - 1, s);
                obj.B1{idx1}.dis = 0;
                obj.B1{idx1}.m = -1;
            end
            
            for k = (obj.windowNum - 1):-1:1
                for s = 0:(obj.clusterNum - 1)
                    idx1 = getNodeIdx(obj, node1VideoIdx, k - 1, s);
                    obj.B1{idx1}.dis = obj.maxDistance;
                    obj.B1{idx1}.m = -1;
                    
                    for m = 0:(obj.clusterNum - 1)
                        idx2 = getNodeIdx(obj, node1VideoIdx, k, m);
                        tmp = obj.B1{idx2}.dis + obj.adjacentMatrix(idx1, idx2);
                        
                        if tmp < obj.B1{idx1}.dis
                            obj.B1{idx1}.dis = tmp;
                            obj.B1{idx1}.m = idx2;
                        end
                    end                  
                end
            end
            
        end
        
        function computeShortestPathForVideo2From0tok(obj)
            k = 0;
            node1VideoIdx = 1;
            
            for s = 0:(obj.clusterNum - 1)
                idx1 = getNodeIdx(obj, node1VideoIdx, k, s);
                obj.A2{idx1}.dis = 0;
                obj.A2{idx1}.m = -1;
            end
            
            for k = 1:(obj.windowNum - 1)
                for s = 0:(obj.clusterNum - 1)
                    idx1 = getNodeIdx(obj, node1VideoIdx, k, s);
                    obj.A2{idx1}.dis = obj.maxDistance;
                    obj.A2{idx1}.m = -1;
                    
                    for m = 0:(obj.clusterNum - 1)
                        idx2 = getNodeIdx(obj, node1VideoIdx, k - 1, m);
                        tmp = obj.A2{idx2}.dis + obj.adjacentMatrix(idx1, idx2);
                        
                        if tmp < obj.A2{idx1}.dis
                            obj.A2{idx1}.dis = tmp;
                            obj.A2{idx1}.m = idx2;
                        end
                    end
                end
            end
        end
        
        function computeShortestPathForVideo2Fromkto0(obj)
            k = obj.windowNum;
            node1VideoIdx = 1;
            
            for s = 0:(obj.clusterNum - 1)
                idx1 = getNodeIdx(obj, node1VideoIdx, k - 1, s);
                obj.B2{idx1}.dis = 0;
                obj.B2{idx1}.m = -1;
            end
            
            for k = (obj.windowNum - 1):-1:1
                for s = 0:(obj.clusterNum - 1)
                    idx1 = getNodeIdx(obj, node1VideoIdx, k - 1, s);
                    obj.B2{idx1}.dis = obj.maxDistance;
                    obj.B2{idx1}.m = -1;
                    
                    for m = 0:(obj.clusterNum - 1)
                        idx2 = getNodeIdx(obj, node1VideoIdx, k, m);
                        tmp = obj.B2{idx2}.dis + obj.adjacentMatrix(idx1, idx2);
                        
                        if tmp < obj.B2{idx1}.dis
                            obj.B2{idx1}.dis = tmp;
                            obj.B2{idx1}.m = idx2;
                        end
                    end                  
                end
            end
        end
        
        function Path = getShortestPathVideo1(obj)
            tmpIdx1 = 0;
            tmpDis = obj.maxDistance;
            node1VideoIdx = 0;
            
            k = obj.windowNum -1;
            for i = 0:(obj.clusterNum - 1)
                idx1 = getNodeIdx(obj, node1VideoIdx, k, i);
                if obj.A1{idx1}.dis < tmpDis
                    tmpDis = obj.A1{idx1}.dis;
                    tmpIdx1 = idx1;
                end
            end
            
            Path.dis = tmpDis;
            Path.nodeIdx = tmpIdx1;
            
            [~, ~, tmpNode1Idx] = getRelativeIdx(obj, tmpIdx1);
            Path.nodeRelative = tmpNode1Idx;
            
            idx1 = tmpIdx1;
            
            for k = obj.windowNum:-1:2
                tmpIdx1 = obj.A1{idx1}.m;
                
                Path.nodeIdx = [tmpIdx1, Path.nodeIdx];
                
                [~, ~, tmpNode1Idx] = getRelativeIdx(obj, tmpIdx1);
                
                Path.nodeRelative = [tmpNode1Idx, Path.nodeRelative];
                
                idx1 = tmpIdx1;
            end
            
            obj.ShortestPath1 = Path;
        end
        
        function Path = getShortestPathVideo1fromkto0(obj)
            tmpIdx1 = 0;
            tmpDis = obj.maxDistance;
            node1VideoIdx = 0;
            
            k = 1;
            for i = 0:(obj.clusterNum - 1)
                idx1 = getNodeIdx(obj, node1VideoIdx, k - 1, i);
                if obj.B1{idx1}.dis < tmpDis
                    tmpDis = obj.B1{idx1}.dis;
                    tmpIdx1 = idx1;
                end
            end
            
            Path.dis = tmpDis;
            Path.nodeIdx = tmpIdx1;
            
            [~, ~, tmpNode1Idx] = getRelativeIdx(obj, tmpIdx1);
            Path.nodeRelative = tmpNode1Idx;
            
            idx1 = tmpIdx1;
            
            for k = 1:obj.windowNum - 1
                tmpIdx1 = obj.B1{idx1}.m;
                
                Path.nodeIdx = [Path.nodeIdx, tmpIdx1];
                
                [~, ~, tmpNode1Idx] = getRelativeIdx(obj, tmpIdx1);
                
                Path.nodeRelative = [Path.nodeRelative, tmpNode1Idx];
                
                idx1 = tmpIdx1;
            end
            
            obj.ShortestPathReverse1 = Path;            
        end
        
        function Path = getShortestPathVideo2(obj)
            tmpIdx1 = 0;
            tmpDis = obj.maxDistance;
            node1VideoIdx = 1;
            
            k = obj.windowNum -1;
            for i = 0:(obj.clusterNum - 1)
                idx1 = getNodeIdx(obj, node1VideoIdx, k, i);
                if obj.A2{idx1}.dis < tmpDis
                    tmpDis = obj.A2{idx1}.dis;
                    tmpIdx1 = idx1;
                end
            end
            
            Path.dis = tmpDis;
            Path.nodeIdx = tmpIdx1;
            
            [~, ~, tmpNode1Idx] = getRelativeIdx(obj, tmpIdx1);
            Path.nodeRelative = tmpNode1Idx;
            
            idx1 = tmpIdx1;
            
            for k = obj.windowNum:-1:2
                tmpIdx1 = obj.A2{idx1}.m;
                
                Path.nodeIdx = [tmpIdx1, Path.nodeIdx];
                
                [~, ~, tmpNode1Idx] = getRelativeIdx(obj, tmpIdx1);
                
                Path.nodeRelative = [tmpNode1Idx, Path.nodeRelative];
                
                idx1 = tmpIdx1;
            end
            
            obj.ShortestPath2 = Path;
        end
        
        function Path = getShortestPathVideo2fromkto0(obj)
            tmpIdx1 = 0;
            tmpDis = obj.maxDistance;
            node1VideoIdx = 1;
            
            k = 1;
            for i = 0:(obj.clusterNum - 1)
                idx1 = getNodeIdx(obj, node1VideoIdx, k - 1, i);
                if obj.B2{idx1}.dis < tmpDis
                    tmpDis = obj.B2{idx1}.dis;
                    tmpIdx1 = idx1;
                end
            end
            
            Path.dis = tmpDis;
            Path.nodeIdx = tmpIdx1;
            
            [~, ~, tmpNode1Idx] = getRelativeIdx(obj, tmpIdx1);
            Path.nodeRelative = tmpNode1Idx;
            
            idx1 = tmpIdx1;
            
            for k = 1:obj.windowNum - 1
                tmpIdx1 = obj.B2{idx1}.m;
                
                Path.nodeIdx = [Path.nodeIdx, tmpIdx1];
                
                [~, ~, tmpNode1Idx] = getRelativeIdx(obj, tmpIdx1);
                
                Path.nodeRelative = [Path.nodeRelative, tmpNode1Idx];
                
                idx1 = tmpIdx1;
            end
            
            obj.ShortestPathReverse2 = Path;
        end
        
        function Path = getShortestPathAcrossSpecificNode(obj, idx)
            % getShortestPathAcrossSpecificNode use absolute index
            [tmpNodeVideoIdx, tmpNodeWindowIdx, tmpNodeIdx] = getRelativeIdx(obj, idx);
            
            if tmpNodeVideoIdx == 0
                A = obj.A1;
                B = obj.B1;
            else
                A = obj.A2;
                B = obj.B2;
            end
            
            Path.dis = A{idx}.dis + B{idx}.dis;
            Path.nodeIdx = idx;
            
            % part 1
            curIdx = idx;            
            Path.nodeIdxRelative = tmpNodeIdx;
            
            for k = tmpNodeWindowIdx:-1:1
                tmpIdx = A{curIdx}.m;
                
                Path.nodeIdx = [tmpIdx, Path.nodeIdx];
                
                [~, ~, tmpNodeIdx] = getRelativeIdx(obj, tmpIdx);
                
                Path.nodeIdxRelative = [tmpNodeIdx, Path.nodeIdxRelative];
                
                curIdx = tmpIdx;
            end
            
            % part 2
            curIdx = idx;
            for k = tmpNodeWindowIdx:obj.windowNum - 2
                tmpIdx = B{curIdx}.m;
                
                Path.nodeIdx = [Path.nodeIdx, tmpIdx];
                
                [~, ~, tmpNodeIdx] = getRelativeIdx(obj, tmpIdx);
                
                Path.nodeIdxRelative = [Path.nodeIdxRelative, tmpNodeIdx];
                
                curIdx = tmpIdx;
            end
        end
        
        function Path = getShortestPathAcrossSpecificNodeRelative(obj, nodeVideoIdx, nodeWinIdx, nodeIdx)
            idx = getNodeIdx(obj, nodeVideoIdx, nodeWinIdx, nodeIdx);
            Path = getShortestPathAcrossSpecificNode(obj, idx);
        end
         
        
        
        
        %% one for all
        function computeAll(obj)
            computeShortestPathForVideo1From0tok(obj);
            computeShortestPathForVideo1Fromkto0(obj);
            
            computeShortestPathForVideo2From0tok(obj);
            computeShortestPathForVideo2Fromkto0(obj);
            
            computeCoupledShortestPathfrom0tok(obj);
            computeCoupledShortestPathfromkto0(obj);
            
            % compute shortest path
            getShortestPath(obj);
            getShortestPathfromkto0(obj);
            
            getShortestPathVideo1(obj);
            getShortestPathVideo1fromkto0(obj);
            
            getShortestPathVideo2(obj);
            getShortestPathVideo2fromkto0(obj);
            
        end
        
    end
    
end

