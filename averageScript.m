inputPrefixes = {'Normal', 'iMAT','GIMME','iMATMachado','GIMMEMachado'};
fileOutputExts = {'FALCON','FBA'};
for i=1:length(inputPrefixes)
    for j=1:length(fileOutputExts)
        load(['analyzeFluxScript' inputPrefixes{i} fileOutputExts{j} '.mat']);
	size(allThoroughStats)
	allThoroughStats(end+1,:) = mean(allThoroughStats,1);
	save(['analyzeFluxScript' inputPrefixes{i} fileOutputExts{j} '.mat'],'allThoroughStats','analyzedCellLines');
    end
end