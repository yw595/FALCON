function overlapSub = getMostCommSub(reducedModel,largeModel)

uniqSubsystems = unique(largeModel.subSystems);
uniqSubsystemsTemp = {};
lengthIntArr = [];
for i=1:length(uniqSubsystems)
    if ~any(strcmp(reducedModel.subSystems,uniqSubsystems{i}))
        uniqSubsystemsTemp{end+1} = uniqSubsystems{i};
        subRxns = largeModel.rxns(strcmp(largeModel.subSystems, uniqSubsystems{i}));
        [~, subRxnIdxs, ~] = intersect(largeModel.rxns,subRxns);
        subMets = largeModel.mets(sum(abs(largeModel.S(:,subRxnIdxs)),2)~=0);
        lengthInt = length(intersect(reducedModel.mets,subMets))/length(subMets);
        lengthIntArr(end+1) = lengthInt;
    end
end

uniqSubsystems = uniqSubsystemsTemp;
[~, sortIdxs] = sort(lengthIntArr);
uniqSubsystems = uniqSubsystems(sortIdxs);

overlapSub = uniqSubsystems{end};
end




        