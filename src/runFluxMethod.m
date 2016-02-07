function fluxes = runFluxMethod(expressionData,expressionIDs,inputName,model,algorithm,expressionSDs)
%initCobraToolbox;
%FBAModel = changeObjective(model,'biomass_reaction');
%FBASoln = optimizeCbModel(model,1);
%v_fba = FBASoln.x;
%run initCobraToolbox to avoid parfor bug
%initCobraToolbox;

if strcmp(algorithm,'EFlux')
    fluxes = call_EFlux(model,expressionIDs, expressionData, model.rxns{strcmp(model.rxnNames,'Biomass')},1);
else        
    % run FALCON
    fluxes = runFALCONStripped(model,expressionIDs,expressionData,expressionSDs,1);
end

end