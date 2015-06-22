%[specificModeliMAT, specificRxnsiMAT] = createTissueSpecificModel(origRecon2,expressionDataHD,1,1,[],'Shlomi');
%[specificModelGIMME, specificRxnsGIMME] = createTissueSpecificModel(origRecon2,expressionDataHD,1,1,[],'GIMME');

% MachadoGIMMEFluxes = call_GIMME(origRecon2, expressionIDsMachado, ...
%     expressionDataMachado, .9, quantile(expressionDataMachado,.25));
% MachadoiMATFluxes = call_iMAT(origRecon2, expressionIDsMachado, expressionDataMachado, ...
%     quantile(expressionDataMachado,.25), quantile(expressionDataMachado,.75), 1);
specificModeliMATMachado = extractSubNetwork(origRecon2,origRecon2.rxns(MachadoiMATFluxes>0));
specificModelGIMMEMachado = extractSubNetwork(origRecon2,origRecon2.rxns(MachadoGIMMEFluxes>0));
save('specificModeliMAT.mat','specificModeliMAT');
save('specificModelGIMME.mat','specificModelGIMME');
save('specificModeliMATMachado.mat','specificModeliMATMachado');
save('specificModelGIMMEMachado.mat','specificModelGIMMEMachado');