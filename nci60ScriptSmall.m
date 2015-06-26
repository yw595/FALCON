inputDirs = {['NCI60Sims' filesep 'nci60mRNA'], ['NCI60Sims' filesep 'nci60prot']};
outputDirs = {'outputSmallNormal', 'outputSmalliMAT', 'outputSmallGIMME', ...
    'outputSmalliMATMachado', 'outputSmallGIMMEMachado'};

load specificModeliMAT.mat;
load specificModelGIMME.mat;
load specificModeliMATMachado.mat;
load specificModelGIMMEMachado.mat;
models = {origRecon2, specificModeliMAT, specificModelGIMME, ...
    specificModeliMATMachado, specificModelGIMMEMachado};

for i=1:length(inputDirs)
    inputDir = inputDirs{i};
    inputFiles = dir(inputDir);

    for j=1:length(outputDirs)
    	outputDir = outputDirs{j};
        if ~exist([inputDir filesep outputDir])
            system(['mkdir ' inputDir filesep outputDir]);
        end

        for k=1:10
            if ~isempty(regexp(inputFiles(k).name, '.csv$'))
                runFALCONStripped(constrainMediumExc(initializeRecon2(models{j})), ...
                [inputDir filesep inputFiles(k).name], 1, [inputDir filesep outputDir]);
            end
        end
    end
end