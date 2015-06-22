inputDirs = {['NCI60Sims' filesep 'nci60mRNA'], ['NCI60Sims' filesep 'nci60prot']};
[cellLinesArray, ~, ~] = readJainTable();

for i=1:length(inputDirs)
    inputDir = inputDirs{i};
    outputDir = 'outputRealSmalliMATMachado';
    if ~exist([inputDir filesep outputDir])
        system(['mkdir ' inputDir filesep outputDir]);
    end

    for k=1:10
    if exist([convertExpressionFileName(cellLinesArray{k}) 'SmallModel.mat'],'file')
       
       load([convertExpressionFileName(cellLinesArray{k}) 'SmallModel.mat'])
       runFALCONStripped(constrainMediumExc(initializeRecon2( specificModeliMATMachado )), [inputDir filesep convertExpressionFileName(cellLinesArray{k}) '.csv'], 1, [inputDir filesep outputDir]);

    end
    end
end