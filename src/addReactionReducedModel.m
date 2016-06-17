function addedModel = addReactionReducedModel(reducedModel,internalIdx,type,biomassMets)

if strcmp(type,'transport')
    reducedModel.rxns{end+1}=['TR_' reducedModel.mets{internalIdx}];
    reducedModel.rev(end+1)=1;
    reducedModel.lb(end+1)=-1000;
    reducedModel.ub(end+1)=1000;
elseif strcmp(type,'demand')
    reducedModel.rxns{end+1}=['DEMAND_' reducedModel.mets{internalIdx}];
    reducedModel.rev(end+1)=0;
    reducedModel.lb(end+1)=0;
    reducedModel.ub(end+1)=1000;
elseif strcmp(type,'biomass')
    reducedModel.rxns{end+1}='BIOMASS';
    reducedModel.rev(end+1)=0;
    reducedModel.lb(end+1)=0;
    reducedModel.ub(end+1)=1000;
end

if strcmp(type,'biomass')
    reducedModel.c(end+1)=1;
    added=0;
    for i=1:length(reducedModel.mets)
        if(sum(strcmp([reducedModel.mets{i}],biomassMets))~=0)
            internalIdx=i;
            if(added==0)
                reducedModel.S(internalIdx,end+1)=-1;
                added=1;
            else
                reducedModel.S(internalIdx,end)=-1;             
            end
        end
    end
else
    reducedModel.c(end+1)=0;
    if strcmp(type,'transport')
        reducedModel.S(internalIdx,end+1)=1;
    else
        reducedModel.S(internalIdx,end+1)=-1;
    end
end

reducedModel.rxnGeneMat(end+1,:)=zeros(1,size(reducedModel.rxnGeneMat,2));
reducedModel.grRules{end+1}='';
reducedModel.subSystems{end+1}='';
reducedModel.rxnNames{end+1}='';
reducedModel.rules{end+1}='';

addedModel = reducedModel;
