function modelToRun = addBiomass(modelRef, modelToRun)
%
% Add biomass pseudoreaction to a model that doesn't have one, based on template model
%
% INPUTS
%       modelRef - template model that has a biomass reaction (assumed origRecon2)
%       modelToRun - model without biomass reaction
%
% OUTPUTS
%       modelToRun - same model as input modelToRun
%       except with added biomass reaction derived from modelRef (see below for details)
%
% Author: Yiping Wang, 2015

% check for reaction that matches name of biomass reaction in origRecon2
if ~any(strcmp(modelToRun.rxns,'biomass_reaction'))

    %find all metabolites that contribute to biomass in modelRef
    biomassMets = modelRef.mets(find(modelRef.S(...
    1:length(modelRef.mets),strcmp( modelRef.rxns,'biomass_reaction' ))~=0));

    % for each of biomassMets, check if modelToRun contains the met
    % if yes, add met to metsList and corresponding stoichiometric coeff to stoichList
    % then add biomass reaction based on metsList and stoichList to modelToRun
    metsList = {}; stoichList = [];
    for l=1:length(biomassMets)
        refIdx = find(strcmp(modelRef.mets,biomassMets{l}));
        toRunIdx = find(strcmp(modelToRun.mets,biomassMets{l}));
        if ~isempty(toRunIdx)
            metsList{end+1} = biomassMets{l};
            stoichList(end+1) = modelRef.S(refIdx,strcmp(modelRef.rxns,'biomass_reaction'));
        end
    end
    modelToRun = addReaction(modelToRun,'biomass_reaction',metsList,stoichList);

end