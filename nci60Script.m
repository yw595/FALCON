inputDirs = {['NCI60Sims' filesep 'nci60mRNA'],['NCI60Sims' filesep 'nci60prot'],['NCI60Sims' filesep 'nci60prot_mRNA']};

for i=1:length(inputDirs)
    inputDir = inputDirs{i}
    outputDir = [inputDir filesep 'output'];
    if ~exist(outputDir)
        system(['mkdir ' outputDir]);
    end

    inputFiles = dir(inputDir);

    for i=1:length(inputFiles)
        if ~isempty(regexp(inputFiles(i).name, '.csv$'))
	    runFALCONStripped(constrainMediumExc(initializeRecon2(origRecon2)), ...
		[inputDir filesep inputFiles(i).name], 1, outputDir);
       end
    end

end