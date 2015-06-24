[cellLinesArray, ~, ~] = readJainTable();

for i=1:length(cellLinesArray)
    if ~exist(['NCI60Sims' filesep 'nci60prot' filesep convertExpressionFileName(cellLinesArray{i}) '.csv'],'file')
	continue;
    end
    totalData = importdata(['NCI60Sims' filesep 'nci60prot' filesep convertExpressionFileName(cellLinesArray{i}) '.csv']);
    expressionDataMachado = totalData.data(2,:);
    expressionIDsMachado = totalData.data(1,:);
    expressionIDsMachado = expressionIDsMachado(~isnan(expressionDataMachado));
    expressionDataMachado = expressionDataMachado(~isnan(expressionDataMachado));

    expressionMedian = median(mean(expressionDataMachado,2));
    expressionData = expressionDataMachado >= expressionMedian;

    makeTissueSpecificModels(cellLinesArray{i}, origRecon2, expressionData, expressionIDsMachado, expressionDataMachado);
end