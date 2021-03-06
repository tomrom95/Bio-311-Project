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

numClusters = 20;

%% Correlation calculations

% Create correlation and cluster tree
corrExpression = corr(expression', 'rows', 'pairwise');
corrDist = 1 - corrExpression;
% clusterTree = linkage(corrDist, 'complete');
% 
% % Show Dendrogram
% figure(1)
% dendrogram(clusterTree,0)
% title('Dendrogram for Hierarchical Clustering');
% drawnow;

%% K Means
% K-means clustering with 16 clusters
[cidx, ctrs] = kmeans(expression, numClusters,...
    'dist','corr',...
    'rep',10,...
    'disp','final');

% Correlation Heat Map
[tmp,ind] = sort(cidx);

corrSquare = corrExpression;
corrSorted = corrSquare(ind, ind);
heatMap = HeatMap(fliplr(corrSorted));
addTitle(heatMap, "Correlation Heat Map");

figure(3)
for c = 1:numClusters
    subplot(numClusters/4,4,c);
    plot(times,expression((cidx == c),:)');
    title(c);
    axis tight
end
drawnow;
suptitle('K-Means Clustering - Expression over time');

% get lists of cluster
cluster1List = getClusterList(cidx, genes, 1);
cluster2List = getClusterList(cidx, genes, 2);
cluster3List = getClusterList(cidx, genes, 3);
cluster4List = getClusterList(cidx, genes, 4);
cluster5List = getClusterList(cidx, genes, 5);
cluster6List = getClusterList(cidx, genes, 6);
cluster7List = getClusterList(cidx, genes, 7);
cluster8List = getClusterList(cidx, genes, 8);
cluster9List = getClusterList(cidx, genes, 9);
cluster10List = getClusterList(cidx, genes, 10);
cluster11List = getClusterList(cidx, genes, 11);
cluster12List = getClusterList(cidx, genes, 12);
cluster13List = getClusterList(cidx, genes, 13);
cluster14List = getClusterList(cidx, genes, 14);
cluster15List = getClusterList(cidx, genes, 15);
cluster16List = getClusterList(cidx, genes, 16);
cluster17List = getClusterList(cidx, genes, 17);
cluster18List = getClusterList(cidx, genes, 18);
cluster19List = getClusterList(cidx, genes, 19);
cluster20List = getClusterList(cidx, genes, 20);

% %% Clustergrams
% % Display clustergram, one for each cluster
% for c = 1:numClusters
%     cgram = clustergram(expression((cidx == c),:),'RowLabels',genes((cidx == c)),...
%     'ColumnLabels',times, 'Cluster', 'column',...
%     'Standardize', 'row');
% 
%     titleLabel = 'Cluster dendrogram and heatmap for cluster %d';
%     titleString = sprintf(titleLabel,c);
%     addTitle(cgram,titleString)
%     drawnow;
%     
% end



