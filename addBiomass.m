function modelToRun = addBiomass(modelRef, modelToRun)

if ~any(strcmp(modelToRun.rxns,'biomass_reaction'))

    biomassMets = modelRef.mets(find(modelRef.S(...
    1:length(modelRef.mets),strcmp( modelRef.rxns,'biomass_reaction' ))~=0));

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