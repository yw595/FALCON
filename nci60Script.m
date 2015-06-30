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
	        totalData = importdata([inputDir filesep cellLinesArray{k} '.csv']);
		expressionIDsMachado = totalData.data(2:end,1);
		expressionDataMachado = totalData.data(2:end,2);
		expressionIDsMachado = expressionIDsMachado(~isnan(expressionDataMachado));
		expressionDataMachado = expressionDataMachado(~isnan(expressionDataMachado));
	        EFluxes = call_EFlux(constrainMediumExc(initializeRecon2(origRecon2)), expressionIDsMachado, expressionDataMachado, 'EX_glc(e)',1);
	    elseif strcmp(outputPrefixes{j},'GXFBA')
	        totalData = importdata([inputDir filesep cellLinesArray{k} '.csv']);
		expressionIDsMachado = totalData.data(2:end,1);
		expressionDataMachado = totalData.data(2:end,2);
		expressionIDsMachado = expressionIDsMachado(~isnan(expressionDataMachado));
		expressionDataMachado = expressionDataMachado(~isnan(expressionDataMachado));
		totalData = importdata([inputDir filesep cellLinesArray{1} '.csv']);
		expressionDataMachadoRef = totalData.data(2:end,2);
		expressionDataMachadoRef = expressionDataMachadoRef(~isnan(expressionDataMachadoRef));
		FBAModel = changeObjective(origRecon2,'biomass_reaction');
		FBASoln = optimizeCbModel(constrainMediumExc(initializeRecon2(FBAModel)),1);
		v_fba = FBASoln.x;
	        GXFBAFluxes = call_GXFBA(constrainMediumExc(initializeRecon2(origRecon2)), expressionIDsMachado, expressionDataMachado, expressionDataMachadoRef, v_fba);
	    else
	    initCobraToolbox;
	    modelToRun = nci60ScriptHelper(origRecon2, outputPrefixes{j}, inputDir, cellLinesArray{k}, INITFilesArray{k}, mCADREFilesArray{k});
	    
            runFALCONStripped(constrainMediumExc(initializeRecon2(modelToRun)), ...
            [inputDir filesep cellLinesArray{k} '.csv'], 1, outputDir);

	    if ~any(strcmp(modelToRun.rxns,'biomass_reaction'))
	        %origRecon2.S(1:length(origRecon2.mets),strcmp( origRecon2.rxns,'biomass_reaction' ))
	        biomassMets = origRecon2.mets(find(origRecon2.S(1:length(origRecon2.mets),strcmp( origRecon2.rxns,'biomass_reaction' ))~=0));
		metsList = {}; stoichList = [];
		for l=1:length(biomassMets)
	            fullIdx = find(strcmp(origRecon2.mets,biomassMets{l}));
	            reducedIdx = find(strcmp(modelToRun.mets,biomassMets{l}));
		    if ~isempty(reducedIdx)
		        metsList{end+1} = biomassMets{l};
			stoichList(end+1) = origRecon2.S(fullIdx,strcmp(origRecon2.rxns,'biomass_reaction'));
		    end
	        end
	        modelToRun = addReaction(modelToRun,'biomass_reaction',metsList,stoichList);
	    end

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