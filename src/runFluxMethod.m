function [fluxes obj cost_rev corr_rho rxn_exp_md_rev] = runFluxMethod(expressionData,expressionIDs, inputName,model,algorithm,expressionSDs,biomassName,overwriteSDs,overwriteMeans)
    
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
elseif strcmp(algorithm,'EFlux')
    if ~exist('biomassName','var')
        biomassName = 'Biomass';
    end
    %disp(expressionIDs)
    %disp(expressionData)
    %disp(model.rxns{strcmp(model.rxnNames,biomassName)})
    fluxes = call_EFlux(model,expressionIDs, expressionData, model.rxns{strcmp(model.rxnNames,biomassName)},1);
    scaleidx = 1;
    while isempty(fluxes)
        fluxes = call_EFlux(model,expressionIDs, expressionData, model.rxns{scaleidx},1);
        scaleidx = scaleidx+1;
    end
    obj = nan;
elseif strcmp(algorithm,'RELATCH')
    modelRELATCH = model;
    expressionIDsRELATCH = expressionIDs;
    for i=1:length(expressionIDs)
        expressionIDsRELATCH{i} = num2str(i);
    end
    for i=1:length(modelRELATCH.genes)
        %modelRELATCH.genes{i} = num2str(find(strcmp(modelRELATCH.genes{i},expressionIDs)));
        modelRELATCH.rules{i} = [ '(x(' num2str(find(strcmp(modelRELATCH.rules{i},expressionIDs))) '))'];    
        %modelRELATCH.grRules{i} = num2str(find(strcmp(modelRELATCH.grRules{i},expressionIDs)));
    end
    fluxes = call_RELATCH(modelRELATCH, expressionIDsRELATCH, expressionData, {}, []);
elseif strcmp(algorithm,'MADE')
    bounds_ref = struct(); bounds_ref.lb = model.lb; bounds_ref.ub = model.ub;
    fluxes = call_MADE(model, expressionIDs, expressionData, expressionData, 0.9, [], bounds_ref);
elseif strcmp(algorithm,'GXFBA')
    fbasol = optimizeCbModel(model);
    fbafluxes = struct();
    fbafluxes.wt_sol = fbasol;
    fbafluxes.wt_minf = model.lb;
    fbafluxes.wt_maxf = model.ub;
    fluxes = call_GXFBA(model, expressionIDs, expressionData, expressionData, fbafluxes);
else      
    % run FALCON
    if ~exist('overwriteSDs','var')
        overwriteSDs = containers.Map;
    end
    if ~exist('overwriteMeans','var')
        overwriteMeans = containers.Map;
    end
    [fluxes falconobj cost_rev corr_rho rxn_exp_md_rev] = runFALCONStripped(model,expressionIDs,expressionData,expressionSDs,1,overwriteSDs,overwriteMeans);
    obj = falconobj;
end

end
