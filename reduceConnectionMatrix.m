function [ reducedMatrix, reducedGenes ] = reduceConnectionMatrix( connectionMatrix, genes )

allZeroIndicies = all(connectionMatrix == 0, 2);
reducedMatrix = connectionMatrix;
reducedMatrix(allZeroIndicies,:) = [];
reducedGenes = genes;
reducedGenes(allZeroIndicies) = [];

end

