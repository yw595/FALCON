inputDirs = {['NCI60Sims' filesep 'nci60prot'], ...
    ['NCI60Sims' filesep 'nci60prot'],['NCI60Sims' filesep 'nci60prot_mRNA']};
outputPrefixes = {'Normal', 'iMAT', 'GIMME', ...
    'iMATMachado', 'GIMMEMachado'};
[cellLinesArray, ~, ~] = readJainTable();

for i=1:length(inputDirs)
    inputDir = inputDirs{i};
    inputFiles = dir(inputDir);
    for j=1:length(outputPrefixes)
        outputDir = [inputDir filesep 'outputPar' outputPrefixes{j}];
        if ~exist(outputDir,'dir')
            system(['mkdir ' outputDir]);
        end
        
        parfor k=1:length(cellLinesArray)
            if strcmp(outputPrefixes{j},'Normal')
                modelToRun = origRecon2;
            else
                load([inputDir filesep 'specificModelsPar' filesep 'specificModel' cellLinesArray{k} outputPrefixes{j} '.mat']);
		eval(['modelToRun = specificModel' outputPrefixes{j} ';']);
            end
            %if ~isempty(regexp(inputFiles(k).name, '.csv$'))
                runFALCONStripped(constrainMediumExc(initializeRecon2(modelToRun)), ...
                    [inputDir filesep cellLinesArray{k} '.csv'], 1, outputDir);
            %end
        end
    end
end