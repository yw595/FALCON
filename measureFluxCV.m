inputDir = ['NCI60Sims' filesep 'nci60mRNA' filesep 'output'];
inputFiles = dir(inputDir);

fluxMatrix = [];
for i=1:length(inputFiles)
    %if ~strcmp(inputFiles(i).name, '.') && ~strcmp(inputFiles(i).name, '..') && ~strcmp(inputFiles(i).name, '.git') && ~strcmp(inputFiles(i).name, '.gitignore')
    if ~isempty(regexp(inputFiles(i).name, '.csv.flux$'))
       fluxData = importdata([inputDir filesep inputFiles(i).name]);
       fluxMatrix(:,end+1) = fluxData.data;
    end
end

fluxMeans = mean(abs(fluxMatrix),2);
%fluxStds = std(fluxMatrix,0,2);
fluxStds = mad(fluxMatrix')';

save('measureFluxCV.mat','fluxMeans','fluxStds');

sortedAbsFluxMatrix = sort(abs(fluxMatrix),'descend');
totalFluxRow = sum(sortedAbsFluxMatrix);
cdfs = [];
rowCount = 1;
while isempty(cdfs) || cdfs(end)<.95
    cdfs(end+1) = mean(sortedAbsFluxMatrix(rowCount,:)./totalFluxRow);
    if length(cdfs)>1
        cdfs(end) = cdfs(end) + cdfs(end-1);
    end
    rowCount = rowCount+1;
end

save('largestFluxContributions.mat','cdfs');