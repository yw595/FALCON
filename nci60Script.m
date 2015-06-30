inputDirs = {['NCI60Sims' filesep 'nci60prot']};%, ...
%    ['NCI60Sims' filesep 'nci60prot'],['NCI60Sims' filesep 'nci60prot_mRNA']};
outputPrefixes = {'Normal', 'iMAT', 'GIMME', ...
    'iMATMachado', 'GIMMEMachado','EFlux','GXFBA'};
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
        
        parfor k=1:length(cellLinesArray)
            if ~strcmp(cellLinesArray{k},'MCF7') && ~strcmp(cellLinesArray{k},'K562')
                if strcmp(outputPrefixes{j},'EFlux')
                    initCobraToolbox;
                    [expressionIDsMachado, expressionDataMachado, ~] = ...
                    readExpressionFile([inputDir filesep cellLinesArray{k} '.csv']);
                    EFluxes = call_EFlux(constrainMediumExc(initializeRecon2(origRecon2)), ...
                    expressionIDsMachado, expressionDataMachado, 'EX_glc(e)',1);
                elseif strcmp(outputPrefixes{j},'GXFBA')
                    initCobraToolbox;
                    [expressionIDsMachado, expressionDataMachado, ~] = ...
                    readExpressionFile([inputDir filesep cellLinesArray{k} '.csv']);
                    [~, expressionDataMachadoRef, ~] = ...
                    readExpressionFile([inputDir filesep cellLinesArray{1} '.csv']);

                    FBAModel = changeObjective(origRecon2,'biomass_reaction');
                    FBASoln = optimizeCbModel(constrainMediumExc(initializeRecon2(FBAModel)),1);
                    v_fba = FBASoln.x;

                    GXFBAFluxes = call_GXFBA(constrainMediumExc(initializeRecon2(origRecon2)), ...
                    expressionIDsMachado, expressionDataMachado, expressionDataMachadoRef, v_fba);
                else
                    initCobraToolbox;
                    modelToRun = nci60ScriptHelper(origRecon2, outputPrefixes{j}, ...
                    inputDir, cellLinesArray{k}, INITFilesArray{k}, mCADREFilesArray{k});
                    
                    runFALCONStripped(constrainMediumExc(initializeRecon2(modelToRun)), ...
                        [inputDir filesep cellLinesArray{k} '.csv'], 1, outputDir);
                    
                    modelToRun = addBiomass(origRecon2,modelToRun);
                    
                    modelToRun = changeObjective(modelToRun,'biomass_reaction');
                    FBASoln = optimizeCbModel(constrainMediumExc(initializeRecon2(modelToRun)),1);
                    v_fba = FBASoln.x;
                    
                    nci60ScriptHelper2([outputDir filesep cellLinesArray{k} '.csv' '_fba_flux.mat'], v_fba);
                    fluxOutputFI = fopen([outputDir filesep cellLinesArray{k} '.csv' '.fba.flux'],'w');
                    for l=1:length(v_fba)
                        fprintf(fluxOutputFI,['R_' modelToRun.rxns{l} ',' num2str(v_fba(l)) '\n']);
                    end
                    fclose(fluxOutputFI);
                end
            end
        end
    end
end