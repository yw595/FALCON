function [fluxes obj] = runFluxMethod(expressionData,expressionIDs, inputName,model,algorithm,expressionSDs,biomassName)
    
if strcmp(algorithm,'FBA')
    %initCobraToolbox;
    if ~exist('biomassName','var')
        biomassName = 'biomass_reaction';
    end
    FBAModel = changeObjective(model,'biomass_reaction');
    FBASoln = optimizeCbModel(model,1);
    v_fba = FBASoln.x;
    fluxes = v_fba
    %run initCobraToolbox to avoid parfor bug
    %initCobraToolbox;
    obj = FBASoln.f;
else if strcmp(algorithm,'EFlux')
    if ~exist('biomassName','var')
        biomassName = 'Biomass';
    end
    fluxes = call_EFlux(model,expressionIDs, expressionData, model.rxns{strcmp(model.rxnNames,'Biomass')},1);
    obj = nan;
else      
    % run FALCON
    [fluxes falconobj] = runFALCONStripped(model,expressionIDs,expressionData,expressionSDs,1);
    obj = falconobj;
end

end