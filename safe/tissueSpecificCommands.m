if ~exist('specificModeliMATHD.mat','file')
    [specificModeliMATHD, specificRxnsiMATHD] = createTissueSpecificModel(origRecon2,expressionDataHD,1,1,[],'Shlomi');
    [specificModelGIMMEHD, specificRxnsGIMMEHD] = createTissueSpecificModel(origRecon2,expressionDataHD,1,1,[],'GIMME');
end

specificFluxesGIMMEHDMachado = call_GIMME(origRecon2, expressionIDsHDMachado, ...
    expressionDataHDMachado, .9, quantile(expressionDataMachado,.25));
specificFluxesiMATHDMachado = call_iMAT(origRecon2, expressionIDsHDMachado, expressionDataHDMachado, ...
    quantile(expressionDataHDMachado,.25), quantile(expressionDataHDMachado,.75), 1);
specificModeliMATHDMachado = extractSubNetwork(origRecon2,origRecon2.rxns(specificFluxesiMATHDMachado>0));
specificModelGIMMEHDMachado = extractSubNetwork(origRecon2,origRecon2.rxns(specificFluxesGIMMEHDMachado>0));

save('specificModeliMATHD.mat','specificModeliMATHD');
save('specificModelGIMMEHD.mat','specificModelGIMMEHD');
save('specificModeliMATHDMachado.mat','specificModeliMATHDMachado');
save('specificModelGIMMEHDMachado.mat','specificModelGIMMEHDMachado');