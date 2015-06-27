inputDirs = {['NCI60Sims' filesep 'nci60prot']};%, ...
%    ['NCI60Sims' filesep 'nci60prot'],['NCI60Sims' filesep 'nci60prot_mRNA']};
outputPrefixes = {'Normal'};%, 'iMAT', 'GIMME', ...
%    'iMATMachado', 'GIMMEMachado'};
[cellLinesArray, ~, ~, ~, ~, ] = readJainTable();
[originTissuesArray INITFilesArray mCADREFilesArray] = makeOriginTissuesArray(cellLinesArray);

for i=1:length(inputDirs)
    inputDir = inputDirs{i};
    inputFiles = dir(inputDir);
    for j=1:length(outputPrefixes)
        outputDir = [inputDir filesep 'output' outputPrefixes{j}];
        if ~exist(outputDir,'dir')
            system(['mkdir ' outputDir]);
        end

        parfor k=1:2%length(cellLinesArray)
	    if strcmp(outputPrefixes{j},'Normal')
                modelToRun = origRecon2;
            else
                load([inputDir filesep 'specificModels' filesep 'specificModel' cellLinesArray{k} outputPrefixes{j} '.mat']);
		eval(['modelToRun = specificModel' outputPrefixes{j} ';']);
            end
            runFALCONStripped(constrainMediumExc(initializeRecon2(modelToRun)), ...
            [inputDir filesep cellLinesArray{k} '.csv'], 1, outputDir);
	    modelToRun = changeObjective(modelToRun,'biomass_reaction');
	    FBASoln = optimizeCbModel(constrainMediumExc(initializeRecon2(modelToRun)),1);
	    v_fba = FBASoln.x;
	    save([outputDir filesep cellLinesArray{k} '.csv' '_fba_flux.mat'], 'v_fba');
	    fluxOutputFI = fopen([outputDir filesep cellLinesArray{k} '.csv' '.flux'],'w');
	    for l=1:length(v_fba)
		fprintf(fluxOutputFI,['R_' modelToRun.rxns{l} ',' num2str(v_fba(l)) '\n']);
	    end
	    fclose(fluxOutputFI);
        end
    end
end