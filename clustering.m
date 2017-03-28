%% Setup
close all
expressionTable = readtable('causton-2001-heat-expression.csv');

% Get expression from table
expression = table2array(expressionTable(:,3:end));
% Get gene names from table
genes = table2array(expressionTable(:,2));

% Filter out any values of '20' (seems to be the NaN value)
nanIndices = any(expression == 20,2);
expression(nanIndices,:) = [];
genes(nanIndices) = [];

colMean = mean(expression(:,1:2),2);
expression = [colMean, expression(:,3:end)];

times = [0,15, 30, 45, 60, 120];

numClusters = 12;

%% Hierarchical clustering
% Hierarchical cluster 
corrDist = pdist(expression, 'corr');
clusterTree = linkage(corrDist, 'complete');
clusters = cluster(clusterTree, 'maxclust', numClusters);

figure(1)
dendrogram(clusterTree,0)
title('Dendrogram for Hierarchical Clustering');
drawnow;

% Plot time series of clusters
figure(2)
for c = 1:numClusters
    subplot(numClusters/4,4,c);
    plot(times,expression((clusters == c),:)');
    axis tight
end
drawnow;
suptitle('Hierarchical Clustering of Profiles');

%% K Means
% K-means clustering with 16 clusters
[cidx, ctrs] = kmeans(expression, numClusters,...
    'dist','corr',...
    'rep',5,...
    'disp','final');
figure(3)
for c = 1:numClusters
    subplot(numClusters/4,4,c);
    plot(times,expression((cidx == c),:)');
    axis tight
end
drawnow;
suptitle('K-Means Clustering of Profiles');

%% Clustergrams
% Display clustergram, one for each cluster
for c = 1:numClusters
    cgram = clustergram(expression((cidx == c),:),'RowLabels',genes((cidx == c)),...
    'ColumnLabels',times, 'Cluster', 'column',...
    'Standardize', 'row');

    titleLabel = 'Cluster dendrogram and heatmap for cluster %d';
    titleString = sprintf(titleLabel,c);
    addTitle(cgram,titleString)
    drawnow;
    
end


