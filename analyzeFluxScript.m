inputPrefixes = {'Normal', 'iMAT','GIMME','iMATMachado','GIMMEMachado'};
[cellLinesArray, ~, ~, ~, ~, ] = readJainTable();
for i=1:length(inputPrefixes)
    inputDir = ['NCI60Sims/nci60prot/output' inputPrefixes{i}];
    inputFiles = dir(inputDir);
    
    fileExts = {'.csv_falcon_flux.mat','.csv_fba_flux.mat'};
    fileOutputExts = {'FALCON','FBA'};
    for k=1:length(fileExts)
        allThoroughStats = [];
        analyzedCellLines = {};
        for j=1:length(cellLinesArray)
            if ~strcmp(cellLinesArray{j},'MCF7') && ~strcmp(cellLinesArray{j},'K562')
                analyzedCellLines{end+1} = cellLinesArray{j};
                if strcmp(inputPrefixes{i},'Normal')
                    modelToRun = origRecon2;
                else
                    load(['NCI60Sims' filesep 'nci60prot' filesep 'specificModelsPar' ...
                    filesep 'specificModel' cellLinesArray{j} inputPrefixes{i} '.mat']);
                    eval(['modelToRun = specificModel' inputPrefixes{i} ';']);
                end
                statsArray = analyzeFlux([inputDir filesep cellLinesArray{j} fileExts{k}], ...
                cellLinesArray{j},modelToRun);
                allThoroughStats(end+1,:) = statsArray(end-3,:);
            end
        end
        allThoroughStats(end+1,:) = mean(allThoroughStats,1);
        
        save(['analyzeFluxScript' inputPrefixes{i} fileOutputExts{k} '.mat'], ...
        'allThoroughStats','analyzedCellLines');
    end
end