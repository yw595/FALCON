function v_Exc = extractExcFlux(model, v_falcon)
%
% Extract exchange fluxes corresponding to CORE data from a flux vector
%
% INPUTS
%       model - model defining what the reaction names are for each flux in v_falcon
%       v_falcon - vector of predicted fluxes
%
% OUTPUTS
%       v_Exc - vector of predicted exchange fluxes, one for each of the CORE metabolites
%
% Author: Yiping Wang, 2015

    [~, jainMetsArray, ~] = readJainTable();
    jainMetsToExcIdxs = loadJainMetsToExcIdxs(jainMetsArray, model);
    v_Exc = [];

    % for each jainMet, either sum up all fluxes in v_falcon 
    % corresponding to it according to jainMetsToExcIdxs,
    % or assign zero if no corresponding fluxes
    for j = 1:length(jainMetsArray)
        jthExcIdxs = jainMetsToExcIdxs(jainMetsArray{j});
	if ~isempty(jthExcIdxs)
            v_Exc(j) = sum(v_falcon(jthExcIdxs));
	else
	    v_Exc(j) = 0;
	end
    end
end