function reducedModel=makeReducedModel(model,subsystemsToAdd,subsystemsToExcludeIdxs,description)

% make manually curated model of 117 reactions that has an FBA solution
curatedRxnNames = makeCuratedRxnNames(model);

%add reactions from 
if ~exist('subsystemsToAdd','var')
    subsystemsToAdd = allSubsystems;
end
if ~exist('subsystemsToExcludeIdxs','var')
    subsystemsToExcludeIdxs = containers.Map;
    for i=1:length(subsystemsToAdd)
        subsystemsToExcludeIdxs(subsystemsToAdd{i}) = [];
    end
end
subsystemRxnNames = makeSubsystemRxnNames(model,subsystemsToAdd,subsystemsToExcludeIdxs);

if (~exist('description','var'))
    description='';
end

% assume subsystemRxnNames is cell array of cell arrays, add rxnName if
% not in curatedRxnNames
reducedRxnNames = curatedRxnNames;
for i=1:length(subsystemRxnNames)
  subsystemRxnNamesSet=subsystemRxnNames{i};
  for j=1:length(subsystemRxnNamesSet)
    if(sum(strcmp(subsystemRxnNamesSet{j},reducedRxnNames))==0)
      reducedRxnNames=[reducedRxnNames; {subsystemRxnNamesSet{j}}];
    end
  end
end
reducedModel = subselectModel(model,description,reducedRxnNames);

% fix four reactions that have weird bounds in origRecon2, NOTE: ALSO PYK LOWER LIM SET TO 1
for i=1:length(reducedModel.rxns)
    if(strcmp(reducedModel.rxns{i},'PYK'))
        reducedModel.lb(i)=1;
	reducedModel.ub(i)=1000;
	reducedModel.rev(i)=0;
    end
    if(strcmp(reducedModel.rxns{i},'PIt2m'))
        reducedModel.lb(i)=-1000;
	reducedModel.ub(i)=1000;
	reducedModel.rev(i)=1;
    end
    if(strcmp(reducedModel.rxns{i},'CYOOm2'))
        reducedModel.lb(i)=-1000;
	reducedModel.ub(i)=1000;
	reducedModel.rev(i)=1;
    end
    if(strcmp(reducedModel.rxns{i},'FTHFDH'))
        reducedModel.lb(i)=-1000;
	reducedModel.ub(i)=1000;
	reducedModel.rev(i)=1;
    end
end

% add transport, drain, and biomass reactions needed for model to function
[transportMets demandMets biomassMets] = makeToAddMets();
for i=1:length(reducedModel.mets)
    internalIdx=i;
    if(sum(strcmp([reducedModel.mets{i}],transportMets))~=0)
        reducedModel = addReactionReducedModel(reducedModel,internalIdx,'transport');
    end
    if(sum(strcmp([reducedModel.mets{i}],demandMets))~=0)
        reducedModel = addReactionReducedModel(reducedModel,internalIdx,'demand');
    end
end
reducedModel = addReactionReducedModel(reducedModel,0,'biomass',biomassMets);

for i=1:length(reducedModel.rxns)
    if(strcmp(reducedModel.rxns{i},'TR_glc_D[c]'))
        reducedModel.lb(i)=-100;
	reducedModel.ub(i)=100;
	reducedModel.rev(i)=1;
    end
end

end