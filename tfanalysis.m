close all
tfTable = readtable('zanton-2004-TFbinding.csv');

genes = table2array(tfTable(:,1));
geneOtherNames = table2array(tfTable(:,2));
for i = 1:length(genes)
   if strcmp(genes(i), '')
      genes(i) = geneOtherNames(i); 
   end
end

tfs = {'TBP', 'TAF1', 'Bdf1', 'Spt3', 'Mot1'};

data = tfTable(:,4:end);
cellTable = table2cell(data);
cellTable = convertCellToDouble(cellTable);
tfMatrix = cell2mat(cellTable);

nanIndices = any(isnan(tfMatrix),2);
tfMatrix(nanIndices,:) = [];
genes(nanIndices) = [];

%% Activators
thresh = 0.6;

connectionMatrix = tfMatrix > thresh;

[activatorMatrix, activatorGenes] = reduceConnectionMatrix(connectionMatrix, genes);

%plotNetwork(connectionMatrix, tfs, activatorGenes);

%count nodes and edges
connectionMatrixT = connectionMatrix';
activatorsPerGene = sum(connectionMatrixT);
numGenesActivated = 0;
for k = 1:length(activatorsPerGene)
    if activatorsPerGene(k) ~= 0
numGenesActivated = numGenesActivated + 1;
    end
end
numEdgesActivation = sum(sum(connectionMatrix));

% Pull outdegrees for each TF and plot histograms for indegrees of each gene
outDegreesActivation = sum(connectionMatrix);
activatorsPerGeneNonzero = activatorsPerGene(activatorsPerGene ~= 0);
figure(1);
histogram(activatorsPerGeneNonzero);
title('Indegrees per gene for activation (genes with no edges excluded)')
xlabel('Indegree (TFs)')
ylabel('Frequency (genes)')
% ***appropriate print command

%% Repressors
thresh = -1.3;

connectionMatrix = tfMatrix < thresh;

[repressorMatrix, repressorGenes] = reduceConnectionMatrix(connectionMatrix, genes);

%plotNetwork(repressorMatrix, tfs, repressorGenes);

%%
%count nodes and edges
connectionMatrixT = connectionMatrix';
repressorsPerGene = sum(connectionMatrixT);
numGenesRepressed = 0;
for k = 1:length(repressorsPerGene)
    if repressorsPerGene(k) ~= 0
numGenesRepressed = numGenesRepressed + 1;
    end
end
numEdgesRepression = sum(sum(connectionMatrix));

% Pull outdegrees for each TF and plot histograms for indegrees of each gene

outDegreesRepression = sum(connectionMatrix);2
repressorsPerGeneNonzero = repressorsPerGene(repressorsPerGene ~= 0);
figure(2);
histogram(repressorsPerGeneNonzero);
title('Indegrees per gene for repression (genes with no edges excluded)')
xlabel('Indegree (TFs)')
ylabel('Frequency (genes)')
% ***appropriate print command


%% Display Linked Genes

disp("Linked Genes");
disp("As activators");
for i = 1:length(tfs)
    disp(tfs(i));
    linkedGeneIndices = activatorMatrix(:,i) == 1;
    linkedGenes = activatorGenes(linkedGeneIndices);
    disp(linkedGenes);
end

% Code to display genes without quotes
% for i = 1:length(linkedGenes)
%     text = sprintf('%s\n',linkedGenes{i});
%     disp(text)
% end

disp("As repressors");
 
for i = 1:length(tfs)
     disp(tfs(i));
     linkedGeneIndices = repressorMatrix(:,i) == 1;
     linkedGenes = repressorGenes(linkedGeneIndices);
     disp(linkedGenes);
end