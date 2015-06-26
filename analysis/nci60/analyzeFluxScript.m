    %inputDir = 'NCI60Sims/nci60prot/outputRealSmalliMATMachado';
    inputDir = 'NCI60Sims/nci60prot/outputSmallNormal';
    inputFiles = dir(inputDir);

    allAverageStats = [];
    analyzedCellLines = {};
    for i=1:length(inputFiles)
    if ~isempty(regexp(inputFiles(i).name, '.csv_falcon_flux.mat'))
        [startIdx endIdx] = regexp(inputFiles(i).name, '.csv_falcon_flux.mat');
        cellLineName = inputFiles(i).name(1:startIdx-1)
	analyzedCellLines{end+1} = cellLineName;
	%load([convertExpressionFileName(cellLineName) 'SmallModel.mat'])
        statsArray = analyzeFlux([inputDir filesep inputFiles(i).name],convertExpressionFileName(cellLineName),origRecon2);%specificModeliMATMachado);
	allAverageStats(end+1,:) = statsArray(end,:);
    end
    end

    save('analyzeFluxScript.mat','allAverageStats','analyzedCellLines');