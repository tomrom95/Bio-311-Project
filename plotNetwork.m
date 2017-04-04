function plotNetwork( connectionMatrix, tfs, genes )
% connectionMatrix: connections with tfs as rows, genes as columns
% tfs: row of tf names
% genes: column of gene names
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


end

