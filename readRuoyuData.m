ruoyuHumanEnsemblToEntrez = importdata('RuoyuEnsemblToEntrez.txt');
ruoyuHumanEnsemblToEntrez.textdata = ruoyuHumanEnsemblToEntrez.textdata(2:end,1);
ruoyuEnsemblMouseToHuman = table2cell(readtable('RuoyuEnsemblMouseToHuman.csv','Format','%s%s'));
if ~exist('ruoyuEnsemblToEntrez','var')
ruoyuEnsemblToEntrez = containers.Map;
for i=1:size(ruoyuEnsemblMouseToHuman,1)
    i
    lineNum = find(strcmp( ruoyuHumanEnsemblToEntrez.textdata(:,1),ruoyuEnsemblMouseToHuman{i,2} ));
    if(length(lineNum)==1 && lineNum > 0)
        if(~isnan(ruoyuHumanEnsemblToEntrez.data(lineNum)))
            ruoyuEnsemblToEntrez(ruoyuEnsemblMouseToHuman{i,1}) = ruoyuHumanEnsemblToEntrez.data(lineNum);
        end
    end
end
end


if ~exist('G1data','var')
    G1data = table2cell(readtable('G1cells.csv'));
end
if ~exist('G2Mdata','var')
    G2Mdata = table2cell(readtable('G2Mcells.csv'));
end
if ~exist('Sdata','var')
    Sdata = table2cell(readtable('Scells.csv'));
end

allEnsemblIDs = importdata('ruoyuMouseEnsemblIDs.txt');

hasEntrezID = zeros(length(allEnsemblIDs),1);
for i=1:length(allEnsemblIDs)
    if isKey(ruoyuEnsemblToEntrez, allEnsemblIDs{i})
        hasEntrezID(i) = 1;
    end
end
hasEntrezID = logical(hasEntrezID);

datasets = {G1data G2Mdata Sdata};
filenames = {'G1FalconData.txt','G2MFalconData.txt','SFalconData.txt'};
for k=1:length(datasets)
    data = datasets{k};
    dataStripped = cell2mat(data(:,2:end));
    dataStripped = dataStripped(:,hasEntrezID);
    ensemblGeneIDArray = allEnsemblIDs(hasEntrezID);
    entrezGeneIDArray = cellfun(@(x) ruoyuEnsemblToEntrez(x), ensemblGeneIDArray);
    meanArray = mean(dataStripped,1); stdArray = std(dataStripped,0,1);
    
    outputFile = fopen(filenames{k},'w');
    fprintf(outputFile,'entrez Gene ID\tmean\tstd\n');
    for i=1:length(entrezGeneIDArray)
        fprintf(outputFile,'%d\t%f\t%f',entrezGeneIDArray(i), meanArray(i), stdArray(i));
        if i~=length(entrezGeneIDArray)
            fprintf(outputFile,'\n');
        end
    end
    fclose(outputFile);
end