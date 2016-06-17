function subModel = subselectModel(model, description, curatedRxnNames, metNames)

% subselect rxns, mets and S
selRxns = ismember(model.rxns,curatedRxnNames);
subS = model.S(:,selRxns);
if ~exist('metNames','var')
    selMets = ~all(subS == 0,2);
else
    selMets = ismember(model.mets,metNames);
end

subS = subS(selMets,:);

subModel.S = subS;
subModel.rxns = model.rxns(selRxns);
subModel.mets = model.mets(selMets);

% subselect b, metNames, Formulas, rev, lb, ub, c, genes (selGenes are any associated with selRxns),
% geneNames, rxnNames, subSystems
if (isfield(model,'b'))
    subModel.b = model.b(selMets);
end
if (isfield(model,'metNames'))
    subModel.metNames = model.metNames(selMets);
end
if (isfield(model,'metFormulas'))
    subModel.metFormulas = model.metFormulas(selMets);
end
if (isfield(model,'metKeggID'))
    subModel.metKeggID = model.metKeggID(selMets);
end
if (isfield(model,'description'))
    subModel.description = description;
end
if (isfield(model,'rev'))
    subModel.rev = model.rev(selRxns);
end
if (isfield(model,'lb'))
    subModel.lb = model.lb(selRxns);
end 
if (isfield(model,'ub'))
    subModel.ub = model.ub(selRxns);
end
if (isfield(model,'c'))
    subModel.c = model.c(selRxns);
end
if (isfield(model,'genes'))
   newRxnGeneMat = model.rxnGeneMat(selRxns,:); 
   selGenes = sum(newRxnGeneMat)' > 0;
   subModel.rxnGeneMat = newRxnGeneMat(:,selGenes);
   subModel.genes = model.genes(selGenes);
   subModel.grRules = model.grRules(selRxns);
   subModel.rules = model.rules(selRxns);
end
if (isfield(model,'geneNames'))
    subModel.geneNameRules = model.geneNameRules(selRxns);
    subModel.geneNames = model.geneNames(selGenes);
end
if (isfield(model,'subSystems'))
    subModel.subSystems = model.subSystems(selRxns);
end
if(isfield(model,'rxnNames'))
    subModel.rxnNames=model.rxnNames(selRxns);
end
if isfield(model,'reversibleModel')
    subModel.reversibleModel=model.reversibleModel;
end