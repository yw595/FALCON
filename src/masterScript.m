model = origRecon2;
paramSets = {};

% subsystemsToAdd = {};
% subsystemsToExcludeIdxs = containers.Map;
% paramSets{end+1,1} = subsystemsToAdd;
% paramSets{end,2} = subsystemsToExcludeIdxs;
% subsystemsToAdd={'Glycolysis/gluconeogenesis'};
% subsystemsToExcludeIdxs = containers.Map;
% subsystemsToExcludeIdxs('Glycolysis/gluconeogenesis') = 34:length(model.rxns(strcmp(model.subSystems,'Glycolysis/gluconeogenesis')));
% paramSets{end+1,1} = subsystemsToAdd;
% paramSets{end,2} = subsystemsToExcludeIdxs;
% subsystemsToAdd={'Glycolysis/gluconeogenesis'};
% subsystemsToExcludeIdxs = containers.Map;
% subsystemsToExcludeIdxs('Glycolysis/gluconeogenesis') = 35:length(model.rxns(strcmp(model.subSystems,'Glycolysis/gluconeogenesis')));
% paramSets{end+1,1} = subsystemsToAdd;
% paramSets{end,2} = subsystemsToExcludeIdxs;
% subsystemsToAdd={'Glycolysis/gluconeogenesis'};
% subsystemsToExcludeIdxs = containers.Map;
% subsystemsToExcludeIdxs('Glycolysis/gluconeogenesis') = [34 39:length(model.rxns(strcmp(model.subSystems,'Glycolysis/gluconeogenesis')))];
% paramSets{end+1,1} = subsystemsToAdd;
% paramSets{end,2} = subsystemsToExcludeIdxs;
% subsystemsToAdd={'Glycolysis/gluconeogenesis'};
% subsystemsToExcludeIdxs = containers.Map;
% subsystemsToExcludeIdxs('Glycolysis/gluconeogenesis') = [34 40:length(model.rxns(strcmp(model.subSystems,'Glycolysis/gluconeogenesis')))];
% paramSets{end+1,1} = subsystemsToAdd;
% paramSets{end,2} = subsystemsToExcludeIdxs;

% dataTable1 = {'xlabels', '1-Curated Model', ...
%               '2-Before Adding UTP Phosphotransferase', ...
%               '3-After Adding UTP Phosphotransferase', ...
%               '4-Before Adding Alternative Hexokinase', ...
%               '5-After Adding Alternative Hexokinase'};
% dataTable2 = {'yvals'};
% dataTable3 = {'Obj'};
% for i=1:length(paramSets)
%     disp(['Param Set ' num2str(i)]);
%     subsystemsToAdd = paramSets{i,1};
%     subsystemsToExcludeIdxs = paramSets{i,2};
%     reducedModel = makeReducedModel(model,subsystemsToAdd,subsystemsToExcludeIdxs);
%     %[falconfluxes falconobj] = runFluxMethod(data,ids,'test',reducedModel,'FALCON',sds,'BIOMASS');
%     [falconfluxes falconobj] = runFluxMethod(data_2,ids_2,'test',reducedModel,'FALCON',stds_2,'BIOMASS');
%     glcLacRat = falconfluxes(strcmp(reducedModel.rxns,'TR_glc_D[c]'))/ ...
%         falconfluxes(strcmp(reducedModel.rxns,'TR_lac_L[c]'));
%     printModel(reducedModel,falconfluxes, ...
%                ['/home/fs01/yw595/reducedFBA' num2str(i) '.txt'], ...
%                ['/home/fs01/yw595/frequentMetsReduced' num2str(i) '.txt']);
%     dataTable2{end+1} = num2str(glcLacRat);
%     dataTable3{end+1} = num2str(falconobj);
% end
% writeData({dataTable1,dataTable2,dataTable3},'/home/fs01/yw595/FALCONDebug.txt','\t');

allSubsystems = {'Glycolysis/gluconeogenesis', 'Pentose phosphate pathway', ...
                                                   'Purine synthesis', 'Pyrimidine synthesis', ...
                                                   'Oxidative phosphorylation','Citric acid cycle', ...
                                                   'Nucleotide interconversion','Folate metabolism'};

% dataTable1 = {'xlabels'};
% dataTable2 = {'yvals'};
% dataTable3 = {'Obj'};

% for i=0:length(allSubsystems)%2:2
%     subsystemsToExcludeIdxs = containers.Map;
%     if i==0
%         subsystemsToAdd = {};
%         dataTable1{end+1} = '1-No Additions';
%     else
%         subsystemsToAdd = allSubsystems(1:i);
%         for j=1:length(subsystemsToAdd)
%             if j==2
%                 subsystemsToExcludeIdxs(subsystemsToAdd{j}) = [];
%                 %subsystemsToExcludeIdxs(subsystemsToAdd{j}) = [39:39];
%             else
%                 subsystemsToExcludeIdxs(subsystemsToAdd{j}) = [];
%             end
%         end
%         dataTable1{end+1} = [num2str(i+1) '-' allSubsystems{i}];
%     end
%     reducedModel = makeReducedModel(model,subsystemsToAdd,subsystemsToExcludeIdxs);
%     %[falconfluxes falconobj] = runFluxMethod(data,ids,'test',reducedModel,'FALCON',sds,'BIOMASS');
%     [falconfluxes falconobj] = runFluxMethod(data_2,ids_2,'test',reducedModel,'FALCON',stds_2,'BIOMASS');
%     glcLacRat = falconfluxes(strcmp(reducedModel.rxns,'TR_glc_D[c]'))/ ...
%         falconfluxes(strcmp(reducedModel.rxns,'TR_lac_L[c]'));
%     %printModel(reducedModel,falconfluxes, ...
%     %           ['/home/fs01/yw595/reducedFBASub' num2str(i) '.txt'], ...
%     %           ['/home/fs01/yw595/frequentMetsReducedSub' num2str(i) '.txt']);
%     printModel(reducedModel,falconfluxes, ...
%                ['/home/fs01/yw595/reducedFBASub' 'Test' '.txt'], ...
%                ['/home/fs01/yw595/frequentMetsReducedSub' 'Test' '.txt']);
%     dataTable2{end+1} = num2str(glcLacRat);
%     dataTable3{end+1} = num2str(falconobj);
% end
% writeData({dataTable1,dataTable2,dataTable3},'/home/fs01/yw595/FALCONSeqGlucose.txt','\t');

if 0
subsystemsToExcludeIdxs = containers.Map;
subsystemsToAdd = {};
reducedModel = makeReducedModel(model,subsystemsToAdd,subsystemsToExcludeIdxs);
[falconfluxes falconobj] = runFluxMethod(data_test,ids_test,'test',reducedModel,'FALCON',stds_test,'BIOMASS');
end

if 0
subsystemsToExcludeIdxs = containers.Map;
subsystemsToAdd = {};
for i=1:length(allSubsystems)
    subsystemsToAdd{end+1} = allSubsystems{i};
    for j=1:length(subsystemsToAdd)
            subsystemsToExcludeIdxs(subsystemsToAdd{j}) = [];
    end
end

reducedModel = makeReducedModel(model,subsystemsToAdd,subsystemsToExcludeIdxs);

reducedModelArr = {reducedModel};
while length(unique(reducedModelArr{end}.subSystems)) ~= length(unique(origRecon2.subSystems))
    mostComm = getMostCommSub(reducedModelArr{end},origRecon2);
    mostComm
    subsystemsToAdd{end+1} = mostComm;
    subsystemsToExcludeIdxs(mostComm) = [];
    reducedModelArr{end+1} = makeReducedModel(model,subsystemsToAdd,subsystemsToExcludeIdxs);
end
end

if 1
ratioArr = [];
parfor i=1:length(reducedModelArr)
    initCobraToolbox;
    [falconfluxes falconobj] = runFluxMethod(data_2,ids_2,'test',reducedModelArr{i},'FALCON',stds_2,'BIOMASS');
    ratioArr(i) = falconfluxes(strcmp(reducedModelArr{i}.rxns,'TR_glc_D[c]'))/ ...
         falconfluxes(strcmp(reducedModelArr{i}.rxns,'TR_lac_L[c]'));
end
end