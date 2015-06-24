function makeTissueSpecificModels(prefix, model, expressionData, expressionIDsMachado, expressionDataMachado)

if ~exist(['smallModels' filesep 'specificModel' prefix 'iMAT.mat'],'file')
    [specificModeliMAT, specificRxnsiMAT] = createTissueSpecificModel(model,expressionData,1,1,[],'Shlomi');
    [specificModelGIMME, specificRxnsGIMME] = createTissueSpecificModel(model,expressionData,1,1,[],'GIMME');
end

specificFluxesGIMMEMachado = call_GIMME(model, expressionIDsMachado, ...
    expressionDataMachado, .9, quantile(expressionDataMachado,.25));
specificFluxesiMATMachado = call_iMAT(model, expressionIDsMachado, expressionDataMachado, ...
    quantile(expressionDataMachado,.25), quantile(expressionDataMachado,.75), 1);
specificModeliMATMachado = extractSubNetwork(model,model.rxns(specificFluxesiMATMachado>0));
specificModelGIMMEMachado = extractSubNetwork(model,model.rxns(specificFluxesGIMMEMachado>0));

save(['smallModels' filesep 'specificModel' prefix 'iMAT.mat'],'specificModeliMAT');
save(['smallModels' filesep 'specificModel' prefix 'GIMME.mat'],'specificModelGIMME');
save(['smallModels' filesep 'specificModel' prefix 'iMATMachado.mat'],'specificModeliMATMachado');
save(['smallModels' filesep 'specificModel' prefix 'GIMMEMachado.mat'],'specificModelGIMMEMachado');

end