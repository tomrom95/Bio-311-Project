close all
tfTable = readtable('zanton-2004-TFbinding.csv');

genes = table2array(tfTable(:,2));
tfs = {'TBP', 'TAF1', 'Bdf1', 'Spt3', 'Mot1'};

data = tfTable(:,4:end);
cellTable = table2cell(data);
cellTable = convertCellToDouble(cellTable);
tfMatrix = cell2mat(cellTable);

nanIndices = any(isnan(tfMatrix),2);
tfMatrix(nanIndices,:) = [];
genes(nanIndices) = [];

thresh = 0.5;

connectionMatrix = tfMatrix > thresh;

allZeroIndicies = all(connectionMatrix == 0, 2);
connectionMatrix(allZeroIndicies,:) = [];
genes(allZeroIndicies) = [];

padding = fliplr(size(connectionMatrix));
symmetric = padarray(connectionMatrix, padding, 'post')';

bg = biograph(symmetric, [tfs genes']);
view(bg);