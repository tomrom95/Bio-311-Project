function [clusterList] = getClusterList(clustInfo,geneList,clustNum)

clusterBoolean = clustInfo == clustNum;
clusterIndices = find(clusterBoolean);
clusterList = geneList(clusterIndices);