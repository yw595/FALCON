% this script creates cell line-specific models using two versions of either the GIMME or iMAT algorithms,
% and outputs .mat files for each specific model to the folder specificModels under each inputDir

% we will look in these inputDirs for expression files ending in .csv for each cell line
inputDirs = {['NCI60Sims' filesep 'nci60prot']};
%    ['NCI60Sims' filesep 'nci60prot'],['NCI60Sims' filesep 'nci60prot_mRNA']};
[cellLinesArray, ~, ~] = readJainTable();

for i=1:length(inputDirs)
    inputDir = inputDirs{i};
    inputFiles = dir(inputDir);

    % make outputDir if necessary
    outputDir = [inputDir filesep 'specificModelsPar' ];
    if ~exist(outputDir,'dir')
        system(['mkdir ' outputDir]);
    end
    
    parfor k=1:length(cellLinesArray)
        % run initCobraToolbox at each iteration, else parfor throws global variable error I don't completely understand yet
        initCobraToolbox;
        cellLinesArray{k}
        
	% read in expression data, then call helper function makeTissueSpecificModels
        [expressionIDsMachado, expressionDataMachado, ~] = ...
        readExpressionFile([inputDir filesep cellLinesArray{k} '.csv']);
        
        makeTissueSpecificModels(outputDir, cellLinesArray{k}, origRecon2, ...
            expressionIDsMachado, expressionDataMachado);
    end
end