readHDData();
load('processedHDData.mat');

expressionDataStriatum = expressionDataHD(:,strcmp(tissues,'striatum'));

makeTissueSpecificModels('HD',origRecon2,expressionDataHD, expressionIDsMachadoHD, expressionDataMachadoHD);