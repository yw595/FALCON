inputDirs = {['NCI60Sims' filesep 'nci60mRNA'],...
    ['NCI60Sims' filesep 'nci60prot'],['NCI60Sims' filesep 'nci60prot_mRNA']};
[cellLinesArray, ~, ~] = readJainTable();

for i=1:length(inputDirs)
    inputDir = inputDirs{i};
    inputFiles = dir(inputDir);
    outputDir = [inputDir filesep 'specificModels' ];
    for k=1:length(cellLinesArray)
        totalData = importdata([inputDir filesep cellLinesArray{k} '.csv']);
        expressionDataMachado = totalData.data(2,:);
        expressionIDsMachado = totalData.data(1,:);
        expressionIDsMachado = expressionIDsMachado(~isnan(expressionDataMachado));
        expressionDataMachado = expressionDataMachado(~isnan(expressionDataMachado));
        
        makeTissueSpecificModels(outputDir, cellLinesArray{k}, origRecon2, ...
            expressionIDsMachado, expressionDataMachado);
    end
end