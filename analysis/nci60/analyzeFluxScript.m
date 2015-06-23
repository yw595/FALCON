    inputDir = 'NCI60Sims/nci60prot/outputRealSmalliMATMachado';
    inputFiles = dir(inputDir);

    for i=1:length(inputFiles)
    if ~isempty(regexp(inputFiles(i).name, '.csv_falcon_flux.mat'))
        [startIdx endIdx] = regexp(inputFiles(i).name, '.csv_falcon_flux.mat');
        cellLineName = inputFiles(i).name(1:startIdx-1)
	load([convertExpressionFileName(cellLineName) 'SmallModel.mat'])
        statsArray = analyzeFlux([inputDir filesep inputFiles(i).name],cellLineName,specificModeliMATMachado);
    end
    end