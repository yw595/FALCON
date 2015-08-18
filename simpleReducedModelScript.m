model = origRecon2;
reducedModel = makeReducedModel(origRecon2);
% subsystems are transferred one at a time from subSystemsToTest
% to Add
subSystemsToTest = {'Glycolysis/gluconeogenesis','Pentose phosphate pathway','Purine synthesis','Pyrimidine synthesis','Citric acid cycle','Oxidative phosphorylation','Exchange/demand reaction','Transport, mitochondrial','Transport, extracellular','Transport, nuclear','Nucleotide interconversion','Fatty acid synthesis','Purine catabolism','Pyrimidine catabolism','Folate metabolism'};
subSystemsToAdd = {};
% by using this block, searching for either glc_D[c] or [e], 
% then using printRxnEq, found five transport rxns
% GLCGLUT2 and GLCt1r, both just exchange, opposite directions
% GLCt2_2, cotransport with 2 protons, GLCt4, cotransport with one sodium
% GLCSGLT1le, cotransport with two sodiums and ATP hydrolysis
glucoseRxns = origRecon2.rxns(origRecon2.S(strcmp(origRecon2.mets,'glc_D[c]'),:)~=0);

glucoseRxnFluxesMatrix = [];
while ~isempty(subSystemsToTest)

    % find largest shared fraction subsystem among subSystemsToTest,
    % (hence ToTest suffix), remove from ToTest and add to ToAdd
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

    % print out largestSharedSubSystem
    largestSharedSubSystem

    % add subSystemsToAdd (hence ToAdd suffix) using makeReducedModel,
    % run FBA with 'BIOMASS' objective, make glucoseRxnFluxes among glucoseRxns present in model
    % add to RxnFluxesMatrix
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

% relevantRxns and RxnFluxesMatrix strip out all nan's, then save
relevantRxns = glucoseRxns(any(~isnan(glucoseRxnFluxesMatrix),1));
glucoseRxnFluxesMatrix = glucoseRxnFluxesMatrix(:,any(~isnan(glucoseRxnFluxesMatrix),1));
save('simpleReducedModelScript.mat','relevantRxns','glucoseRxnFluxesMatrix','subSystemsToAdd');

% ???, possibly just print out earlier
reducedModel = makeReducedModel(origRecon2);
reducedModel = changeObjective(reducedModel,'BIOMASS',1);
reducedModelFBA = optimizeCbModel(reducedModel);
for i=1:length(reducedModelFBA.x)
  %disp(sprintf('%s\t%f',reducedModel.rxns{i},reducedModelFBA.x(i)))
end
%reducedModelFBA.x(strcmp(reducedModel.rxns,'HEX1'))