inputDirs = {['NCI60Sims' filesep 'nci60mRNA'],...
    ['NCI60Sims' filesep 'nci60prot'],['NCI60Sims' filesep 'nci60prot_mRNA']};
outputPrefixes = {'Normal', 'iMAT', 'GIMME', ...
    'iMATMachado', 'GIMMEMachado'};
[cellLinesArray, ~, ~] = readJainTable();

for i=1:length(inputDirs)
    inputDir = inputDirs{i};
    inputFiles = dir(inputDir);
    for j=1:length(outputPrefixes)
        outputDir = [inputDir filesep 'output' outputPrefixes{j}];
        if ~exist(outputDir,'dir')
            system(['mkdir ' outputDir]);
        end
        
        for k=1:length(cellLinesArray)
            if strcmp(outputPrefixes{j},'Normal')
                modelToRun = origRecon2;
            else
                load([inputDir filesep cellLinesArray{k} outputPrefixes{j} 'Model.mat']);
            end
            if ~isempty(regexp(inputFiles(k).name, '.csv$'))
                runFALCONStripped(constrainMediumExc(initializeRecon2(modelToRun)), ...
                    [inputDir filesep inputFiles(k).name], 1, outputDir);
            end
        end
    end
end