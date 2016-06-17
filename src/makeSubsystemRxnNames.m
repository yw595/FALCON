function subsystemRxnNames = makeSubsystemRxnNames(model,subsystemsToAdd,subsystemsToExcludeIdxs)

subsystemNamesToRxns = containers.Map;
subsystemRxnNames = {};
for i=1:length(subsystemsToAdd)
    subsystemNamesToAdd = model.rxns(strcmp(model.subSystems,subsystemsToAdd{i}));
    subsystemNamesToAdd = subsystemNamesToAdd(setdiff(1:length(subsystemNamesToAdd),subsystemsToExcludeIdxs(subsystemsToAdd{i})));
    subsystemNamesToRxns(subsystemsToAdd{i}) = subsystemNamesToAdd;
    subsystemRxnNames{end+1} = subsystemNamesToRxns(subsystemsToAdd{i});
end

% glycolysisRxnNames=model.rxns(strcmp(model.subSystems,'Glycolysis/gluconeogenesis'));
% pentosePhosphateRxnNames=model.rxns(strcmp(model.subSystems,'Pentose phosphate pathway'));
% purineSynthesisRxnNames=model.rxns(strcmp(model.subSystems,'Purine synthesis'));
% pyrimidineSynthesisRxnNames=model.rxns(strcmp(model.subSystems,'Pyrimidine synthesis'));
% citricAcidCycleRxnNames=model.rxns(strcmp(model.subSystems,'Citric acid cycle'));
% oxidativePhosphorylationRxnNames=model.rxns(strcmp(model.subSystems,'Oxidative phosphorylation'));
% exchangeDemandRxnNames=model.rxns(strcmp(model.subSystems,'Exchange/demand reaction'));
% transportMitoRxnNames=model.rxns(strcmp(model.subSystems,'Transport, mitochondrial'));
% transportExRxnNames=model.rxns(strcmp(model.subSystems,'Transport, extracellular'));
% transportNuclearRxnNames=model.rxns(strcmp(model.subSystems,'Transport, nuclear'));
% nucleotideInterconversionRxnNames=model.rxns(strcmp(model.subSystems,'Nucleotide interconversion'));
% fattyAcidSynthesisRxnNames=model.rxns(strcmp(model.subSystems,'Fatty acid synthesis'));
% purineCatabolismRxnNames=model.rxns(strcmp(model.subSystems,'Purine catabolism'));
% pyrimidineCatabolismRxnNames=model.rxns(strcmp(model.subSystems,'Pyrimidine catabolism'));
% folateMetabolismRxnNames=model.rxns(strcmp(model.subSystems,'Folate metabolism'));

%excludeRxnNames={exchangeDemandRxnNames; transportExRxnNames; purineCatabolismRxnNames; pyrimidineCatabolismRxnNames};
%glycolysisRxnNames(34)
%glycolysisRxnNames(39)
%glycolysisRxnNames = glycolysisRxnNames([1:33]);
%toAddRxnNames={glycolysisRxnNames};%; pentosePhosphateRxnNames; purineSynthesisRxnNames; pyrimidineSynthesisRxnNames; oxidativePhosphorylationRxnNames; citricAcidCycleRxnNames; nucleotideInterconversionRxnNames; folateMetabolismRxnNames};