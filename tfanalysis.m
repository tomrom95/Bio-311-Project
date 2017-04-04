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
padX = padding(1);
padY = padding(2);
symmetric = padarray(connectionMatrix, [padX 0], 'pre');
symmetric = padarray(symmetric, [0 padY], 'post');

tfGraph = biograph(symmetric', [ tfs genes' ], 'LayoutType', 'radial');
view(tfGraph);