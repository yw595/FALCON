[cellLinesArray, ~, ~] = readJainTable();

for i=1:length(cellLinesArray)
    if ~exist(['NCI60Sims' filesep 'nci60prot' filesep convertExpressionFileName(cellLinesArray{i}) '.csv'],'file')
	continue;
    end
    totalData = importdata(['NCI60Sims' filesep 'nci60prot' filesep convertExpressionFileName(cellLinesArray{i}) '.csv']);
    expressionDataMachado = totalData.data(2,:);
    expressionIDsMachado = totalData.data(1,:);
    expressionIDsMachado = expressionIDsMachado(~isnan(expressionDataMachado));
    expressionDataMachado = expressionDataMachado(~isnan(expressionDataMachado));
    
    MachadoGIMMEFluxes = call_GIMME(origRecon2, expressionIDsMachado, ...
    expressionDataMachado, .9, quantile(expressionDataMachado,.25));
    MachadoiMATFluxes = call_iMAT(origRecon2, expressionIDsMachado, expressionDataMachado, ...
    quantile(expressionDataMachado,.25), quantile(expressionDataMachado,.75), 1);

    specificModeliMATMachado = extractSubNetwork(origRecon2,origRecon2.rxns(MachadoiMATFluxes>0));
    specificModelGIMMEMachado = extractSubNetwork(origRecon2,origRecon2.rxns(MachadoGIMMEFluxes>0));

    save([convertExpressionFileName(cellLinesArray{i}) 'SmallModel.mat'],'specificModeliMATMachado','specificModelGIMMEMachado','MachadoiMATFluxes','MachadoGIMMEFluxes');
end