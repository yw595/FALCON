model = origRecon2;
reducedModel = makeReducedModel(origRecon2);
% subsystems are transferred one at a time from subSystemsToTest
% to Add
subSystemsToTest = {'Glycolysis/gluconeogenesis','Pentose phosphate pathway','Purine synthesis','Pyrimidine synthesis','Citric acid cycle','Oxidative phosphorylation','Exchange/demand reaction','Transport, mitochondrial','Transport, extracellular','Transport, nuclear','Nucleotide interconversion','Fatty acid synthesis','Purine catabolism','Pyrimidine catabolism','Folate metabolism'};
%subSystemsToTest = [subSystemsToTest, setdiff(unique(origRecon2.subSystems), subSystemsToTest)];
subSystemsToTest2 = setdiff(unique(origRecon2.subSystems), subSystemsToTest);
disp(subSystemsToTest)
subSystemsToAdd = {};
% by using this block, searching for either glc_D[c] or [e], 
% then using printRxnEq, found five transport rxns
% GLCGLUT2 and GLCt1r, both just exchange, opposite directions
% GLCt2_2, cotransport with 2 protons, GLCt4, cotransport with one sodium
% GLCSGLT1le, cotransport with two sodiums and ATP hydrolysis
glucoseRxns = origRecon2.rxns(origRecon2.S(strcmp(origRecon2.mets,'glc_D[c]'),:)~=0);
% add these two here so they show up in bar graph
glucoseRxns{end+1} = 'EX_glc(e)'; glucoseRxns{end+1} = 'TR_glc_D[c]';
consumingGlucoseRxns = origRecon2.rxns(origRecon2.S(strcmp(origRecon2.mets,'glc_D[c]'),:)<0);
producingGlucoseRxns = origRecon2.rxns(origRecon2.S(strcmp(origRecon2.mets,'glc_D[c]'),:)>0);

glucoseRxnFluxesMatrix = [];
while ~isempty(subSystemsToTest) || ~isempty(subSystemsToTest2)

    if isempty(subSystemsToTest)
        subSystemsToTest = subSystemsToTest2;
	subSystemsToTest2 = {};
    end
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
    % following arrays sum up exchange, consuming and producing fluxes
    glucoseExchangeArray(end+1) = reducedModelFBA.x(strcmp(reducedModel.rxns,'TR_glc_D[c]'));
    if any(strcmp(reducedModel.rxns,'EX_glc(e)'))
        glucoseExchangeArray(end) = glucoseExchangeArray(end) + reducedModelFBA.x(strcmp(reducedModel.rxns,'EX_glc(e)'));
    end
    [~, consumingIdxs, ~] = intersect(consumingGlucoseRxns, glucoseRxns);
    glucoseConsumingArray(end+1) = sum(glucoseRxnFluxes(consumingIdxs));
    [~, producingIdxs, ~] = intersect(producingGlucoseRxns, glucoseRxns);
    if ~isempty(glucoseRxnFluxes(producingIdxs))
        glucoseProducingArray(end+1) = sum(glucoseRxnFluxes(producingIdxs));
    else
        glucoseProducingArray(end+1) = nan;
    end
end

% relevantRxns and RxnFluxesMatrix strip out all nan's, then save
relevantRxns = glucoseRxns(any(~isnan(glucoseRxnFluxesMatrix),1));
glucoseRxnFluxesMatrix = glucoseRxnFluxesMatrix(:,any(~isnan(glucoseRxnFluxesMatrix),1));
save('simpleReducedModelScript.mat','relevantRxns','glucoseRxnFluxesMatrix','glucoseExchangeArray','glucoseConsumingArray','glucoseProducingArray','subSystemsToAdd');

% ???, possibly just print out earlier
reducedModel = makeReducedModel(origRecon2);
reducedModel = changeObjective(reducedModel,'BIOMASS',1);
reducedModelFBA = optimizeCbModel(reducedModel);
for i=1:length(reducedModelFBA.x)
  %disp(sprintf('%s\t%f',reducedModel.rxns{i},reducedModelFBA.x(i)))
end
%reducedModelFBA.x(strcmp(reducedModel.rxns,'HEX1'))