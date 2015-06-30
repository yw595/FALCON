function modelToRun = addBiomass(modelRef, modelToRun)

if ~any(strcmp(modelToRun.rxns,'biomass_reaction'))
    biomassMets = modelRef.mets(find(modelRef.S(...
    1:length(modelRef.mets),strcmp( modelRef.rxns,'biomass_reaction' ))~=0));
    metsList = {}; stoichList = [];
    for l=1:length(biomassMets)
        fullIdx = find(strcmp(modelRef.mets,biomassMets{l}));
        reducedIdx = find(strcmp(modelToRun.mets,biomassMets{l}));
        if ~isempty(reducedIdx)
            metsList{end+1} = biomassMets{l};
            stoichList(end+1) = modelRef.S(fullIdx,strcmp(modelRef.rxns,'biomass_reaction'));
        end
    end
    modelToRun = addReaction(modelToRun,'biomass_reaction',metsList,stoichList);
end

end