function modelToRun = nci60ScriptHelper(origRecon2, outputPrefix, inputDir, cellLine, INITFile, mCADREFile)
    if strcmp(outputPrefix,'Normal')
	modelToRun = origRecon2;
    elseif any(strcmp(outputPrefix,{'iMAT','GIMME','iMATMachado','GIMMEMachado'}))
	load([inputDir filesep 'specificModelsPar' filesep 'specificModel' cellLine outputPrefix '.mat']);
	eval(['modelToRun = specificModel' outputPrefix ';']);
    elseif strcmp(outputPrefix,'INIT')
	if ~strcmp(INITFile,'na')
	    modelToRun = readCbModel(['DownloadedTissueSpecific' filesep 'INITFiles' filesep INITFile '.xml']);
	end
	load(['DownloadedTissueSpecific' filesep 'INITFiles' filesep 'HMRdatabase2_00.mat']);
	HMRModel = ihuman;
	[~, intersectRxnIdxs, ~] = union(HMRModel.rxns, cellfun(@(x) ['HMR_' x], modelToRun.rxns));
	modelToRun.grRules = HMRModel.grRules(intersectRxnIdxs);
	modelToRun.rxnGeneMat = HMRModel.rxnGeneMat(intersectRxnIdxs,:);
	modelToRun.genes = HMRModel.genes;
	ensemblEntrezData = importdata('AllHumanEnsemblToEntrez.txt');
	ensemblEntrezData.textdata = ensemblEntrezData(2:end,:);
	[~, intersectRxnIdxsA, intersectRxnIdxsB] = union(modelToRun.genes, ensemblEntrezData.textdata);
	modelToRun.genes{intersectRxnIdxsA} = num2str(ensemblEntrezData.data(intersectRxnIdxsB));
    else
	if ~strcmp(mCADREFile,'na')
	    load(['DownloadedTissueSpecific' filesep 'mCADREFiles' filesep mCADREFile '.mat']);
	    eval(['modelToRun = ' mCADREFile ';']);
	end
    end

    if strcmp(outputPrefix, 'INIT') || strcmp(outputPrefix,'mCADRE')
        [~, intersectMetIdxsA, intersectMetIdxsB] = union(origRecon2.metFormulas,modelToRun.metFormulas);
	modelToRun.mets{intersectMetIdxsA} = num2str(origRecon2.mets{intersectMetIdxsB});
	modelToRun.metNames{intersectMetIdxsA} = num2str(origRecon2.metNames{intersectMetIdxsB});
    end
end