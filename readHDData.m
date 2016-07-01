haveProcessedData = 0;
if exist('HDDataProcessed.mat','file')
    haveProcessedData = 1;
end

if ~haveProcessedData
    [junk1, junk2, rawMetadata] = xlsread( ...
    '6_month_allelic_series_data_mRNA_with_metadata2014_21.xlsx','mRNA Metadata');
    mouseIDs = rawMetadata(2:end,1);
    tissues = rawMetadata(2:end,3);
    QLengths = rawMetadata(2:end,4);
    sexes = rawMetadata(2:end,6);

    [junk1, ~, rawData] = xlsread( ...
    '6_month_allelic_series_data_mRNA_with_metadata2014_21.xlsx','mRNA FPKM data');
    geneIDs = rawData(2:end,1);
    expressionDataHD = cell2mat(rawData(2:end,2:end));
    save('HDDataProcessed.mat','mouseIDs','tissues','QLengths','sexes','geneIDs','expressionDataHD');
end