% this script examines the simulated flux predictions for selected output folders for each algorithm in NCI60Sims/nci60prot
% it runs analyzeFlux on each one, stores the statistics when looking at all CORE fluxes
% averages these across all cell lines in an algorithm, and outputs them in a set of .mat files

% we will look at folders of the form NCI60Sims/nci60prot/output(inputPrefixes{i})
inputPrefixes = {'EFlux','GXFBA'};%{'Normal', 'iMAT','GIMME','iMATMachado','GIMMEMachado','EFlux','GXFBA'};
[cellLinesArray, ~, ~, ~, ~, ] = readJainTable();
for i=1:length(inputPrefixes)
    inputDir = ['NCI60Sims/nci60prot/output' inputPrefixes{i}];
    inputFiles = dir(inputDir);
    
    % for EFlux and GXFBA, their .mat files end with .csv_eflux_flux.mat or 
    % .csv_gxfba_flux.mat, and we will use EFLUX or GXFBA in their output filenames
    % for all other inputPrefixes, we have both FALCON and FBA simulations to consider
    if strcmp(inputPrefixes{i},'EFlux')
        fileExts = {'.csv_eflux_flux.mat'};
	fileOutputExts = {'EFLUX'};
    elseif strcmp(inputPrefixes{i},'GXFBA')
        fileExts = {'.csv_gxfba_flux.mat'};
	fileOutputExts = {'GXFBA'};
    else
        fileExts = {'.csv_falcon_flux.mat','.csv_fba_flux.mat'};
	fileOutputExts = {'FALCON','FBA'};
    end

    for k=1:length(fileExts)
        allThoroughStats = [];
        analyzedCellLines = {};
        parfor j=1:length(cellLinesArray)
	    % skip two cell lines where GIMME and iMAT failed to make models
            if ~strcmp(cellLinesArray{j},'MCF7') && ~strcmp(cellLinesArray{j},'K562')
	        disp(cellLinesArray{j})
                analyzedCellLines{j} = cellLinesArray{j};

		% load either origRecon2, or a specific model using helper function, depending on inputPrefix
                if strcmp(inputPrefixes{i},'Normal') || strcmp(inputPrefixes{i},'EFlux') || strcmp(inputPrefixes{i},'GXFBA')
                    modelToRun = origRecon2;
                else
                    modelToRun = analyzeFluxScriptHelper(cellLinesArray{j}, inputPrefixes{i});
                end

		% run analyzeFlux, save the fourth row from the end of statsArray,
		% which contains results when looking at all CORE fluxes
                statsArray = analyzeFlux([inputDir filesep cellLinesArray{j} fileExts{k}], ...
                cellLinesArray{j},modelToRun);
                allThoroughStats(j,:) = statsArray(end-3,:);
            end
        end

	% append row for stats averaged across all cell lines
	% save this average, and individual stats for each cell line, 
	%in .mat file, one for each algorithm we look at
        allThoroughStats(end+1,:) = mean(allThoroughStats,1);
        save(['analyzeFluxScript' inputPrefixes{i} fileOutputExts{k} '.mat'], ...
        'allThoroughStats','analyzedCellLines');
    end
end