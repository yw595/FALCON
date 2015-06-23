haveReadXLS = 1;
if exist('rawMetadata','var') && exist('rawData','var')
    haveReadXLS = 0;
end

if(haveReadXLS)
    [junk1, junk2, rawMetadata] = xlsread( ...
    '6_month_allelic_series_data_mRNA_with_metadata2014_21.xlsx','mRNA Metadata');
    mouseIDs = rawMetadata(2:end,1);
    tissues = rawMetadata(2:end,3);
    QLengths = rawMetadata(2:end,4);
    sexes = rawMetadata(2:end,6);

    [junk1, ~, rawData] = xlsread( ...
    '6_month_allelic_series_data_mRNA_with_metadata2014_21.xlsx','mRNA FPKM data');
    geneIDs = rawData(2:end,1);
end
expressionData = cell2mat(rawData(2:end,2:end));
expressionData = expressionData(:,strcmp(tissues,'striatum'));

ensemblMouseToHumanData = textscan(fopen('All Mouse Ensembl To Human.csv'),'%s%s','Delimiter',',','HeaderLines',1);
ensemblMouseToHuman = containers.Map(ensemblMouseToHumanData{1}, ...
    ensemblMouseToHumanData{2});
ensemblHumanToEntrezData = importdata('All Human Ensembl To Entrez.txt');
ensemblHumanToEntrez = containers.Map(ensemblHumanToEntrezData.textdata(2:end,1), ...
    ensemblHumanToEntrezData.data(:,1));
ensemblMouseToEntrez = containers.Map;
ensemblMouseToHumanKeys = keys(ensemblMouseToHuman);
ensemblHumanToEntrezKeys = keys(ensemblHumanToEntrez);
for i=1:length(ensemblMouseToHumanKeys)
    if isKey(ensemblMouseToHuman, ensemblMouseToHumanKeys{i})
        ensemblHumanID = ensemblMouseToHuman(ensemblMouseToHumanKeys{i});
        if isKey(ensemblHumanToEntrez, ensemblHumanID)
            entrezID = ensemblHumanToEntrez(ensemblHumanID);
            if ~isempty(entrezID)
                ensemblMouseToEntrez(ensemblMouseToHumanKeys{i}) = entrezID;
            end
        end
    end
end

ensemblMouseToEntrezKeys = keys(ensemblMouseToEntrez);
[intersection,intersectIDsA,~] = intersect(geneIDs,ensemblMouseToEntrezKeys);
expressionData = expressionData(intersectIDsA,:);
expressionMedian = median(mean(expressionData,2));
%expressionMedian = median(reshape(expressionData,numel(expressionData),1));
%presentAbsentVector = mean(expressionData,2) >= expressionMedian;

expressionDataHD.Locus = [];
expressionDataHD.Data = [];
for i=1:length(intersection)
    if ~isnan(ensemblMouseToEntrez(intersection{i}))
        expressionDataHD.Locus(end+1) = ensemblMouseToEntrez(intersection{i});
        expressionDataHD.Data(end+1) = mean(expressionData(i,:),2);
    end
end
expressionDataMachado = expressionDataHD.Data;
expressionDataHD.Data = expressionDataHD.Data >= expressionMedian;
expressionIDsMachado = expressionDataHD.Locus;