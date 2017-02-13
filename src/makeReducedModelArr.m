function [reducedModelArr subsystemsAddedArr] = makeReducedModelArr(reducedModel,fullModel,useOrigRxns)

if ~exist('useOrigRxns','var')
    useOrigRxns = 0;
end

reducedModelArr = {reducedModel};
subsystemsAddedArr = {'None'};
subsystemsToExcludeIdxs = containers.Map;
subsystemsToAdd = {};
if exist('setSubsystems','var')
    for i=1:length(setSubsystems)
        subsystemsToAdd{end+1} = setSubsystems{i};
        for j=1:length(subsystemsToAdd)
            subsystemsToExcludeIdxs(subsystemsToAdd{j}) = [];
        end
        if useOrigRxns
            reducedModelArr{end+1} = makeReducedModel(fullModel,subsystemsToAdd,subsystemsToExcludeIdxs,'',reducedModel.rxns);
        else
            reducedModelArr{end+1} = makeReducedModel(fullModel,subsystemsToAdd,subsystemsToExcludeIdxs);
        end
        reducedModelArr{end+1} = reducedModel;
        subsystemsAddedArr{end+1} = setSubsystems{i};
    end
end

while length(unique(reducedModelArr{end}.subSystems)) ~= length(unique(fullModel.subSystems))
    mostComm = getMostCommSub(reducedModelArr{end},fullModel);
    mostComm
    
    subsystemsToAdd{end+1} = mostComm;
    subsystemsToExcludeIdxs(mostComm) = [];
    if useOrigRxns
        reducedModelArr{end+1} = makeReducedModel(fullModel,subsystemsToAdd,subsystemsToExcludeIdxs,'',reducedModel.rxns);
    else
        reducedModelArr{end+1} = makeReducedModel(fullModel,subsystemsToAdd,subsystemsToExcludeIdxs);
    end
    subsystemsAddedArr{end+1} = mostComm;    
end

for i=1:length(subsystemsAddedArr)
    subsystemsAddedArr{i} = strrep(subsystemsAddedArr{i},'/','_');
    subsystemsAddedArr{i} = strrep(subsystemsAddedArr{i},' ','_');
    subsystemsAddedArr{i} = strrep(subsystemsAddedArr{i},',','_');
end

end