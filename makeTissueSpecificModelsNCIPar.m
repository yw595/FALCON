inputDirs = {['NCI60Sims' filesep 'nci60prot']};
%    ['NCI60Sims' filesep 'nci60prot'],['NCI60Sims' filesep 'nci60prot_mRNA']};
[cellLinesArray, ~, ~] = readJainTable();

for i=1:length(inputDirs)
    inputDir = inputDirs{i};
    inputFiles = dir(inputDir);
    outputDir = [inputDir filesep 'specificModelsPar' ];
    if ~exist(outputDir,'dir')
        system(['mkdir ' outputDir]);
    end

    %indices = [21,22,23,24,39,40,41];
    indices = [24,40,41];
    parfor l=1:3%1:length(cellLinesArray)
        %changeCobraSolver('gurobi5','MILP');
        k = indices(l);
        initCobraToolbox;
        cellLinesArray{k}
	k
	if strcmp(cellLinesArray{k}, 'IGROV1')
	    disp('HEREEEEEEEEE')
	end
        totalData = importdata([inputDir filesep cellLinesArray{k} '.csv']);
        expressionIDsMachado = totalData.data(2:end,1);
        expressionDataMachado = totalData.data(2:end,2);
        expressionIDsMachado = expressionIDsMachado(~isnan(expressionDataMachado));
        expressionDataMachado = expressionDataMachado(~isnan(expressionDataMachado));
	
	size(expressionDataMachado)
	size(expressionIDsMachado)
        makeTissueSpecificModels(outputDir, cellLinesArray{k}, origRecon2, ...
            expressionIDsMachado, expressionDataMachado);
    end
end