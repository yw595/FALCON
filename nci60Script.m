% this script runs various flux prediction algorithms (defined in outputPrefixes)
% on the expression files located in each inputDir (currently just NCI60Sims/nci60prot),
% and outputs both .mat and .csv.flux files into the appropriate output dirs

% we will look in these inputDirs for expression files ending in .csv for each cell line
inputDirs = {['NCI60Sims' filesep 'nci60prot']};%, ...
%    ['NCI60Sims' filesep 'nci60prot'],['NCI60Sims' filesep 'nci60prot_mRNA']};
% we will output to folders of the form NCI60Sims/nci60prot/output(outputPrefixes{i})
outputPrefixes = {'RELATCH'};%{'Normal', 'iMAT', 'GIMME', ...
%    'iMATMachado', 'GIMMEMachado','EFlux','GXFBA'};
[cellLinesArray, ~, ~, ~, ~, ] = readJainTable();
[originTissuesArray INITFilesArray mCADREFilesArray] = makeOriginTissuesArray(cellLinesArray);

for i=1:length(inputDirs)
    inputDir = inputDirs{i};
    inputFiles = dir(inputDir);
    for j=1:length(outputPrefixes)
        outputDir = [inputDir filesep 'output' outputPrefixes{j}];
	%make outputDir if necessary
        if ~exist(outputDir,'dir')
            system(['mkdir ' outputDir]);
        end
        
        for k=1:length(cellLinesArray)
	    % skip two cell lines where GIMME and iMAT failed to make models
            if ~strcmp(cellLinesArray{k},'MCF7') && ~strcmp(cellLinesArray{k},'K562')

	        % in what follows, run initCobraToolbox at beginning of each iteraction,
	        % otherwise parfor throws an error involving global variables I don't completely understand yet
	        initCobraToolbox;
		% read in common expressiondata and ids for everything except final else clause
		[expressionIDsMachado, expressionDataMachado, ~] = ...
                    readExpressionFile([inputDir filesep cellLinesArray{k} '.csv']);
		%reference FBA soln for everything except final else clause
		FBAModel = changeObjective(origRecon2,'biomass_reaction');
                FBASoln = optimizeCbModel(constrainMediumExc(initializeRecon2(FBAModel)),1);
                v_fba = FBASoln.x;
	        if strcmp(outputPrefixes{j},'RELATCH')
		    tempIDs = {};
		    for z =1:length(expressionIDsMachado)
		        tempIDs{z} = num2str(expressionIDsMachado(z));
		    end
                    RELATCHFluxes = call_RELATCH(constrainMediumExc(initializeRecon2(origRecon2)), ...
                    tempIDs, expressionDataMachado, {}, []);
		    nci60ScriptHelper2([outputDir filesep cellLinesArray{k} '.csv' '_relatch_flux.mat'], [outputDir filesep cellLinesArray{k} '.csv' '.relatch.flux'], RELATCHFluxes, origRecon2);
		elseif strcmp(outputPrefixes{j},'MADE')
		    [expressionIDsMachadoRef, expressionDataMachadoRef, ~] = ...
                    readExpressionFile([inputDir filesep cellLinesArray{1} '.csv']);
		    [intersectIDs, intersectIdxsA, intersectIdxsB] = intersect(expressionIDsMachado, expressionIDsMachadoRef);
		    intersectDataA = expressionDataMachado(intersectIdxsA);
		    intersectDataB = expressionDataMachadoRef(intersectIdxsB);
		    bounds_ref = struct();
		    bounds_ref.lb = origRecon2.lb;
		    bounds_ref.ub = origRecon2.ub;
		    MADEFluxes = call_MADE(constrainMediumExc(initializeRecon2(origRecon2)), ...
		    num2str(intersectIDs), intersectDataA, intersectDataB, 0.9, bounds_ref);
                elseif strcmp(outputPrefixes{j},'EFlux')
                    EFluxes = call_EFlux(constrainMediumExc(initializeRecon2(origRecon2)), ...
                    expressionIDsMachado, expressionDataMachado, 'EX_glc(e)',1);
		    nci60ScriptHelper2([outputDir filesep cellLinesArray{k} '.csv' '_eflux_flux.mat'], [outputDir filesep cellLinesArray{k} '.csv' '.eflux.flux'], EFluxes, origRecon2);
                elseif strcmp(outputPrefixes{j},'GXFBA')
                    [expressionIDsMachadoRef, expressionDataMachadoRef, ~] = ...
                    readExpressionFile([inputDir filesep cellLinesArray{1} '.csv']);
		    [intersectIDs, intersectIdxsA, intersectIdxsB] = intersect(expressionIDsMachado, expressionIDsMachadoRef);
		    intersectDataA = expressionDataMachado(intersectIdxsA);
		    intersectDataB = expressionDataMachadoRef(intersectIdxsB);

                    GXFBAFluxes = call_GXFBA(constrainMediumExc(initializeRecon2(origRecon2)), ...
                    num2str(intersectIDs), intersectDataA, intersectDataB, GXFBA);
		    nci60ScriptHelper2([outputDir filesep cellLinesArray{k} '.csv' '_gxfba_flux.mat'], [outputDir filesep cellLinesArray{k} '.csv' '.gxfba.flux'], GXFBAFluxes, origRecon2);
                else

		    % load appropriate model depending on inputPrefix
                    modelToRun = nci60ScriptHelper(origRecon2, outputPrefixes{j}, ...
                    inputDir, cellLinesArray{k}, INITFilesArray{k}, mCADREFilesArray{k});
                    
		    % run FALCON
                    runFALCONStripped(constrainMediumExc(initializeRecon2(modelToRun)), ...
                        [inputDir filesep cellLinesArray{k} '.csv'], 1, outputDir);
                    
		    % add biomass reaction if necessary, then change objective to that reaction,
		    % and run FBA
                    modelToRun = addBiomass(origRecon2,modelToRun);
                    modelToRun = changeObjective(modelToRun,'biomass_reaction');
                    FBASoln = optimizeCbModel(constrainMediumExc(initializeRecon2(modelToRun)),1);
                    v_fba = FBASoln.x;
                    
		    % save .mat and .csv.flux files for fba using helper function
                    nci60ScriptHelper2([outputDir filesep cellLinesArray{k} '.csv' '_fba_flux.mat'], [outputDir filesep cellLinesArray{k} '.csv' '.fba.flux'], v_fba);
                end
            end
        end
    end
end