function newModel = prepModel(oldModel)

testRxnSubset = 0; isMouse = 0; isYeast = 0; isHuman = 0;
if(strcmp(oldModel.description,'Mus Musculus iMM1415 model'))
    isMouse = 1;
elseif(strcmp(oldModel.description,'LeeYeast.xml'))
    isYeast = 1;
elseif(strcmp(oldModel.description,'origRecon2'))
    isHuman = 1;
end
newModel = oldModel;

if isMouse
    newModel.c = zeros(length(newModel.rxns),1);
end

if isMouse
for i=1:length(newModel.grRules)
    oldRule = newModel.grRules{i};
    if isMouse
        if length(newModel.grRules{i})==0
            newModel.grRules{i} = '';
        end
    end
    [startIdx endIdx] = regexp(oldRule,'\(\S+\)');
    newRule = '';
    if length(startIdx)>0
        newRule = [newRule oldRule(1:startIdx(1)-1)];
        for j=1:length(startIdx)-1
            newRule = [newRule oldRule(startIdx(j)+1:endIdx(j)-1)];
            newRule = [newRule oldRule(endIdx(j)+1:startIdx(j+1)-1)];
        end
        newRule = [newRule oldRule(startIdx(end)+1:endIdx(end)-1)];
        newRule = [newRule oldRule(endIdx(end)+1:end)];
    end
    if ~isempty(strfind(newRule,'and')) || ~isempty(strfind(newRule,'or'))
        newRule = ['(' newRule ')'];
    end
    newModel.grRules{i} = newRule;
end
end

if testRxnSubset
    cutoff1 = 1;
    cutoff2 = 10;
    newModel.S = newModel.S(:,cutoff1:cutoff2);
    newModel.rxns = newModel.rxns(cutoff1:cutoff2);
    newModel.rev = newModel.rev(cutoff1:cutoff2);
    newModel.lb = newModel.lb(cutoff1:cutoff2);
    newModel.ub = newModel.ub(cutoff1:cutoff2);
    newModel.subSystems = newModel.subSystems(cutoff1:cutoff2);
    newModel.rxnGeneMat = newModel.rxnGeneMat(cutoff1:cutoff2,:);
    newModel.grRules = newModel.grRules(cutoff1:cutoff2);
    newModel.rules = newModel.rules(cutoff1:cutoff2);
    newModel.c = newModel.c(cutoff1:cutoff2);
end

if isMouse
    if ~isempty(find(strcmp(newModel.rxns,'PDE4')))
        newModel.grRules{strcmp(newModel.rxns,'PDE4')} = '';
    end
    if ~isempty(find(strcmp(newModel.rxns,'PDE1')))
        newModel.grRules{strcmp(newModel.rxns,'PDE1')} = '';
    end
    if ~isempty(find(strcmp(newModel.rxns,'PIt7ir')))
        newModel.grRules{strcmp(newModel.rxns,'PIt7ir')} = '';
    end
    if ~isempty(find(strcmp(newModel.rxns,'UAGDP')))
        newModel.grRules{strcmp(newModel.rxns,'UAGDP')} = '';
    end
end

end