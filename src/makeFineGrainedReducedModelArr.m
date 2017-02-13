function [reducedModelArr2,reactionsAddedArr2,ratioArr2,distArr2,corrArr2,avgCostArr2,costArr2,expArr2,subsystemsToAdd2] = makeFineGrainedReducedModelArr(baseSubIdx,randomize,model,reducedModelArr,ratioArr,distArr,corrArr,avgCostArr,costArr,expArr,subsystemsToAdd)

if 1 % make finegrained reducedModelArr2
reducedModelArr2 = {reducedModelArr{baseSubIdx}};
reactionsAddedArr2 = {'None'};
if exist('ratioArr','var')
    ratioArr2 = [ratioArr(baseSubIdx)];
end
if exist('distArr','var')
    distArr2 = {distArr{baseSubIdx}};
    corrArr2 = [corrArr(baseSubIdx)];
    avgCostArr2 = [avgCostArr(baseSubIdx)];
    costArr2 = {costArr{baseSubIdx}};
    expArr2 = {expArr{baseSubIdx}};
    subsystemsToAdd2 = subsystemsToAdd(1:baseSubIdx);
end

%subsystemsToExcludeIdxs2 = containers.Map;
%for i=1:length(subsystemsToAdd2)-1
%    subsystemsToExcludeIdxs2(subsystemsToAdd2{i}) = [];
%end

reactionsToAdd = reducedModelArr{baseSubIdx+1}.rxnNames(strcmp(reducedModelArr{baseSubIdx+1}.subSystems,subsystemsToAdd{baseSubIdx}));
if randomize
    reactionsToAdd = reactionsToAdd(randperm(length(reactionsToAdd)));
end
%excludeArr = {};
subsystemsToExcludeArr = {};
for i=0:length(reactionsToAdd)
    tempSub = containers.Map;
    for j=1:length(subsystemsToAdd2)-1
        tempSub(subsystemsToAdd2{j}) = [];
    end
    %realI = length(reactionsToAdd)-i;
    %tempSub(subsystemsToAdd2{end}) = 1:realI;
    tempSub(subsystemsToAdd2{end}) = i+1:length(reactionsToAdd);
    subsystemsToExcludeArr{i+1} = tempSub;
    %tempSub = subsystemsToExcludeArr{i+1};
    if i~= 0
        reactionsAddedArr2{i+1} = reactionsToAdd{i};
    end
end

end

if 1 % does actual work of making reducedModelArr2
parfor i=1:length(reactionsToAdd)
    reducedModel2 = makeReducedModel(model,subsystemsToAdd2, subsystemsToExcludeArr{i+1});
    reducedModelArr2{i+1} = reducedModel2;
end
end

for i=1:length(reactionsAddedArr2)
    reactionsAddedArr2{i} = strrep(reactionsAddedArr2{i},'/','_');
    reactionsAddedArr2{i} = strrep(reactionsAddedArr2{i},' ','_');
    reactionsAddedArr2{i} = strrep(reactionsAddedArr2{i},',','_');
    reactionsAddedArr2{i} = strrep(reactionsAddedArr2{i},'''','_');
end