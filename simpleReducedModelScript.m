model = origRecon2;
reducedModel = makeReducedModel(origRecon2);
subSystemsToTest = {'Glycolysis/gluconeogenesis','Pentose phosphate pathway','Purine synthesis','Pyrimidine synthesis','Citric acid cycle','Oxidative phosphorylation','Exchange/demand reaction','Transport, mitochondrial','Transport, extracellular','Transport, nuclear','Nucleotide interconversion','Fatty acid synthesis','Purine catabolism','Pyrimidine catabolism','Folate metabolism'};
subSystemsToAdd = {};
glucoseRxns = origRecon2.rxns(origRecon2.S(strcmp(origRecon2.mets,'glc_D[c]'),:)~=0);

glucoseRxnFluxesMatrix = [];
while ~isempty(subSystemsToTest)
largestSharedFraction = 0;
largestSharedSubSystem = '';
for i=1:length(subSystemsToTest)
    subSystemMets = model.mets(sum(abs(model.S(:,strcmp(model.subSystems,subSystemsToTest{i}))),2)~=0);
    sharedFraction = sum(ismember( reducedModel.mets,subSystemMets ))/length(subSystemMets);
    if sharedFraction > largestSharedFraction
        largestSharedFraction = sharedFraction;
	largestSharedSubSystem = subSystemsToTest{i};
    end
end
subSystemsToTest(strcmp( subSystemsToTest,largestSharedSubSystem )) = [];
subSystemsToAdd{end+1} = largestSharedSubSystem;

largestSharedSubSystem

rxnsToAdd = {};
for i=1:length(subSystemsToAdd)
    rxnsToAdd{end+1} = origRecon2.rxns(strcmp( origRecon2.subSystems,subSystemsToAdd(i) ));
end
reducedModel = makeReducedModel(origRecon2,rxnsToAdd);
reducedModel = changeObjective(reducedModel,'BIOMASS',1);
reducedModelFBA = optimizeCbModel(reducedModel);
glucoseRxnFluxes = [];
for i=1:length(glucoseRxns)
    if any(strcmp(reducedModel.rxns,glucoseRxns{i}))
        glucoseRxnFluxes(i) = reducedModelFBA.x(strcmp(reducedModel.rxns,glucoseRxns{i}));
    else
        glucoseRxnFluxes(i) = nan;
    end
end
glucoseRxnFluxesMatrix(end+1,:) = glucoseRxnFluxes;

end

relevantRxns = glucoseRxns(any(~isnan(glucoseRxnFluxesMatrix),1));
glucoseRxnFluxesMatrix = glucoseRxnFluxesMatrix(:,any(~isnan(glucoseRxnFluxesMatrix),1));
save('simpleReducedModelScript.mat','relevantRxns','glucoseRxnFluxesMatrix','subSystemsToAdd');

reducedModel = makeReducedModel(origRecon2);
reducedModel = changeObjective(reducedModel,'BIOMASS',1);
reducedModelFBA = optimizeCbModel(reducedModel);
for i=1:length(reducedModelFBA.x)
  %disp(sprintf('%s\t%f',reducedModel.rxns{i},reducedModelFBA.x(i)))
end
%reducedModelFBA.x(strcmp(reducedModel.rxns,'HEX1'))