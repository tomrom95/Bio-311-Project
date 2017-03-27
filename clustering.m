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

%% Hierarchical clustering
% Hierarchical luster with 16 clusters
corrDist = pdist(expression, 'corr');
clusterTree = linkage(corrDist, 'average');
clusters = cluster(clusterTree, 'maxclust', 16);

% Plot time series of clusters
figure
for c = 1:16
    subplot(4,4,c);
    plot(times,expression((clusters == c),:)');
    axis tight
end
suptitle('Hierarchical Clustering of Profiles');

%% K Means
% K-means clustering with 16 clusters
[cidx, ctrs] = kmeans(expression, 16,...
    'dist','corr',...
    'rep',5,...
    'disp','final');
figure
for c = 1:16
    subplot(4,4,c);
    plot(times,expression((cidx == c),:)');
    axis tight
end
suptitle('K-Means Clustering of Profiles');

% Display clustergram
figure
clustergram(expression,'RowLabels',genes,...
    'ColumnLabels',times, 'Cluster', 'column',...
    'Standardize', 'row')

