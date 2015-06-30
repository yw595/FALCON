inputPrefixes = {'Normal', 'iMAT','GIMME','iMATMachado','GIMMEMachado'};
for i=1:length(inputPrefixes)
    inputDir = ['NCI60Sims/nci60prot/output' inputPrefixes{i}];
    inputFiles = dir(inputDir);
    
    fileExts = {'.csv_falcon_flux.mat','.csv_fba_flux.mat'};
    fileOutputExts = {'FALCON','FBA'};
    for k=1:length(fileExts)
        allThoroughStats = [];
	analyzedCellLines = {};
	for j=1:length(inputFiles)
	    inputFiles(j)
	    %if strcmp(inputFiles(j).name,'UACC_257.csv_falcon_flux.mat') || strcmp(inputFiles(j).name,'OVCAR_8.csv_falcon_flux.mat') || strcmp(inputFiles(j).name,'UACC_257.csv_fba_flux.mat') || strcmp(inputFiles(j).name,'OVCAR_8.csv_fba_flux.mat')
	    if ~isempty(regexp(inputFiles(j).name, fileExts{k}))
		[startIdx endIdx] = regexp(inputFiles(j).name, fileExts{k});
		cellLineName = inputFiles(j).name(1:startIdx-1)
		analyzedCellLines{end+1} = cellLineName;
		if strcmp(inputPrefixes{i},'Normal')
		    modelToRun = origRecon2;
		else
		    load(['NCI60Sims' filesep 'nci60prot' filesep 'specificModelsPar' filesep 'specificModel' cellLineName inputPrefixes{i} '.mat']);
		    eval(['modelToRun = specificModel' inputPrefixes{i} ';']);
		end
		statsArray = analyzeFlux([inputDir filesep inputFiles(j).name],convertExpressionFileName(cellLineName),modelToRun);
		%statsArray
		allThoroughStats(end+1,:) = statsArray(end-3,:);
	    end
	    %end
	end

	save(['analyzeFluxScript' inputPrefixes{i} fileOutputExts{k} '.mat'],'allThoroughStats','analyzedCellLines');
    end
end