function makeTissueSpecificModels(outputDir, prefix, model, expressionIDsMachado, expressionDataMachado)

expressionMedian = median(mean(expressionDataMachado,2));
expressionData.Data = expressionDataMachado >= expressionMedian;
expressionData.Locus = expressionIDsMachado;
            
if ~exist([outputDir filesep 'specificModel' prefix 'iMAT.mat'],'file')
    [specificModeliMAT, specificRxnsiMAT] = createTissueSpecificModel(model,expressionData,1,1,[],'Shlomi');
    [specificModelGIMME, specificRxnsGIMME] = createTissueSpecificModel(model,expressionData,1,1,[],'GIMME');
    save([outputDir filesep 'specificModel' prefix 'iMAT.mat'],'specificModeliMAT');
    save([outputDir filesep 'specificModel' prefix 'GIMME.mat'],'specificModelGIMME');
end

specificFluxesGIMMEMachado = call_GIMME(model, expressionIDsMachado, ...
    expressionDataMachado, .9, quantile(expressionDataMachado,.25));
specificFluxesiMATMachado = call_iMAT(model, expressionIDsMachado, expressionDataMachado, ...
    quantile(expressionDataMachado,.25), quantile(expressionDataMachado,.75), 1);
specificModeliMATMachado = extractSubNetwork(model,model.rxns(specificFluxesiMATMachado>0));
specificModelGIMMEMachado = extractSubNetwork(model,model.rxns(specificFluxesGIMMEMachado>0));

save([outputDir filesep 'specificModel' prefix 'iMATMachado.mat'],'specificModeliMATMachado');
save([outputDir filesep 'specificModel' prefix 'GIMMEMachado.mat'],'specificModelGIMMEMachado');

end