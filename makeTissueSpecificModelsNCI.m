inputDirs = {['NCI60Sims' filesep 'nci60mRNA'], ...
    ['NCI60Sims' filesep 'nci60prot'],['NCI60Sims' filesep 'nci60prot_mRNA']};
[cellLinesArray, ~, ~] = readJainTable();

for i=1:length(inputDirs)
    inputDir = inputDirs{i};
    inputFiles = dir(inputDir);
    outputDir = [inputDir filesep 'specificModels' ];
    if ~exist(outputDir,'dir')
        system(['mkdir ' outputDir]);
    end
    for k=1:length(cellLinesArray)
        totalData = importdata([inputDir filesep cellLinesArray{k} '.csv']);
        expressionIDsMachado = totalData.data(2:end,1);
        expressionDataMachado = totalData.data(2:end,2);
        expressionIDsMachado = expressionIDsMachado(~isnan(expressionDataMachado));
        expressionDataMachado = expressionDataMachado(~isnan(expressionDataMachado));
	
        makeTissueSpecificModels(outputDir, cellLinesArray{k}, origRecon2, ...
            expressionIDsMachado, expressionDataMachado);
    end
end