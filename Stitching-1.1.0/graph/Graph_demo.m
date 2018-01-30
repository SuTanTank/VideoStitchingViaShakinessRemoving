clear;
clc;
close all;

dbstop error;

graph = Graph(4, 3);
% graph = Graph(2, 2); % windowNum = 24, clusterNum = 5 

%% set weights
load('share/graph.mat');

% graph.generateRNDWeights();

% video1
% graph.setEdgeWeight(1,5,5);
% graph.setEdgeWeight(2,5,3);
% 
% graph.setEdgeWeight(5,7,6);
% graph.setEdgeWeight(5,8,7);
% 
% 
% graph.setEdgeWeight(7,11,8);
% graph.setEdgeWeight(8,12,6);
% 
% % video2
% graph.setEdgeWeight(15,17,1);
% graph.setEdgeWeight(17,21,5);
% graph.setEdgeWeight(21,23,7);
% 
% % video 1 and video 2
% graph.setEdgeWeight(1,15,2);
% graph.setEdgeWeight(2,15,10);
% 
% graph.setEdgeWeight(5,17,2);
% 
% graph.setEdgeWeight(8,21,2);
% 
% graph.setEdgeWeight(12,23,3);

%%
% graph.computeCoupledShortestPathfrom0tok();
% graph.computeCoupledShortestPathfromkto0();
% graph.computeShortestPathForVideo1From0tok();
% graph.computeShortestPathForVideo2From0tok();
% graph.computeShortestPathForVideo1Fromkto0();
% graph.computeShortestPathForVideo2Fromkto0();

graph.computeAll();

Path = graph.getShortestPath();
PathReverse = graph.getShortestPathfromkto0();

PathEnforceSingleVideo = graph.getShortestPathAcrossSpecificNode(24);
PathEnforceSingleVideoRelative = graph.getShortestPathAcrossSpecificNodeRelative(1, 1, 0);


% call computeCoupledShortestPathfrom0tok() and
% computeCoupledShortestPathfromkto0()first before enforce two nodes
% PathEnforce = graph.getShortestPathAcrossSpecificNodes(1, 123);
% PathEnforceRelative = graph.getShortestPathAcrossSpecificNodesRelative(0, 2, 1, 1, 2, 2);

% this function is used for evaluation, do not use it for large scale 
% problems, as the running time increase exponentially.
% PathBF = graph.getShortestPathBF(); 

% [videoIdx, winIdx, nodeIdx] = graph.getRelativeIdx(7);