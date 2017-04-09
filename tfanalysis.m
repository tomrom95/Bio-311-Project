close all
tfTable = readtable('zanton-2004-TFbinding.csv');

genes = table2array(tfTable(:,1));
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

%% Repressors
thresh = -1.3;

connectionMatrix = tfMatrix < thresh;

[repressorMatrix, repressorGenes] = reduceConnectionMatrix(connectionMatrix, genes);

%plotNetwork(repressorMatrix, tfs, repressorGenes);

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

%% Display Linked Genes

disp("Linked Genes");
disp("As activators");
for i = 1:length(tfs)
    disp(tfs(i));
    linkedGeneIndices = activatorMatrix(:,i) == 1;
    linkedGenes = activatorGenes(linkedGeneIndices);
    disp(linkedGenes);
end

disp("As repressors");

for i = 1:length(tfs)
    disp(tfs(i));
    linkedGeneIndices = repressorMatrix(:,i) == 1;
    linkedGenes = repressorGenes(linkedGeneIndices);
    disp(linkedGenes);
end