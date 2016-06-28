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
[falconfluxes falconobj] = runFluxMethod(data_257,ids_257,'test',reducedModel,'FALCON',stds_257,'BIOMASS');
%[falconfluxes falconobj] = runFluxMethod(data_test,ids_test,'test',reducedModel,'FALCON',stds_test,'BIOMASS');
end

if 0
reducedModelArr = {reducedModel};
subsystemsAddedArr = {''};
subsystemsToExcludeIdxs = containers.Map;
subsystemsToAdd = {};
for i=1:length(allSubsystems)
    subsystemsToAdd{end+1} = allSubsystems{i};
    for j=1:length(subsystemsToAdd)
            subsystemsToExcludeIdxs(subsystemsToAdd{j}) = [];
    end
    reducedModel = makeReducedModel(model,subsystemsToAdd,subsystemsToExcludeIdxs);
    reducedModelArr{end+1} = reducedModel;
    subsystemsAddedArr{end+1} = allSubsystems{i};
end

%reducedModelArr = {reducedModel};
while length(unique(reducedModelArr{end}.subSystems)) ~= length(unique(origRecon2.subSystems))
    mostComm = getMostCommSub(reducedModelArr{end},origRecon2);%, 'Exchange/demand reaction');
    mostComm
    %if ~strcmp(mostComm,'Exchange/demand reaction')
    subsystemsToAdd{end+1} = mostComm;
    subsystemsToExcludeIdxs(mostComm) = [];
    reducedModelArr{end+1} = makeReducedModel(model,subsystemsToAdd,subsystemsToExcludeIdxs);
    subsystemsAddedArr{end+1} = mostComm;
    %end
end
end

%[falconfluxes falconobj] = runFluxMethod(data_257,ids_257,'test',reducedModel,'FALCON',stds_257,'BIOMASS');

for i=1:length(subsystemsAddedArr)
    subsystemsAddedArr{i} = strrep(subsystemsAddedArr{i},'/','_');
    subsystemsAddedArr{i} = strrep(subsystemsAddedArr{i},' ','_');
    subsystemsAddedArr{i} = strrep(subsystemsAddedArr{i},',','_');
end

if 0
ratioArr = [];
distArr = {};
parfor i=1:length(reducedModelArr)
    initCobraToolbox;
    [falconfluxes falconobj] = runFluxMethod(data_257,ids_257,'test',reducedModelArr{i},'FALCON',stds_257,'BIOMASS');
    ratioArr(i) = falconfluxes(strcmp(reducedModelArr{i}.rxns,'TR_glc_D[c]'))/ ...
         -falconfluxes(strcmp(reducedModelArr{i}.rxns,'TR_lac_L[c]'));
    if any(strcmp(reducedModelArr{i}.rxns,'EX_glc(e)')) && any(strcmp(reducedModelArr{i}.rxns,'EX_lac_L(e)'))
        ratioArr(i) = ratioArr(i) - falconfluxes(strcmp(reducedModelArr{i}.rxns,'EX_glc(e)'))/ ...
         falconfluxes(strcmp(reducedModelArr{i}.rxns,'EX_lac_L(e)'));
    end
    distArr{i} = falconfluxes;
end
end

if 0
parfor i=1:length(reducedModelArr)
    i
    printModel(reducedModelArr{i},distArr{i}, ...
               ['/home/fs01/yw595/reducedFBASubAll' subsystemsAddedArr{i} num2str(i) '.txt'], ...
               ['/home/fs01/yw595/frequentMetsReducedSubAll' subsystemsAddedArr{i} num2str(i) '.txt']);
end
end

if 1
reducedModelArr2 = {reducedModelArr{76}};
reactionsAddedArr2 = {''};
ratioArr2 = [ratioArr(76)];
distArr2 = {distArr{76}};
subsystemsToAdd2 = subsystemsToAdd(1:76);
subsystemsToExcludeIdxs2 = containers.Map;
for i=1:length(subsystemsToAdd2)-1
    subsystemsToExcludeIdxs2(subsystemsToAdd2{i}) = [];
end
reactionsToAdd = reducedModelArr{77}.rxnNames(strcmp(reducedModelArr{77}.subSystems,subsystemsToAdd{76}));
excludeArr = {};
subsystemsToExcludeArr = {};
for i=0:length(reactionsToAdd)
    tempSub = containers.Map;
    for j=1:length(subsystemsToAdd2)-1
        tempSub(subsystemsToAdd2{j}) = [];
    end
    realI = length(reactionsToAdd)-i;
    tempSub(subsystemsToAdd2{end}) = 1:realI;
    subsystemsToExcludeArr{i+1} = tempSub;
    tempSub = subsystemsToExcludeArr{i+1};
end

if 1
parfor i=0:length(reactionsToAdd)
    i
    realI = length(reactionsToAdd)-i;
    subsystemsToExcludeIdxs3 = subsystemsToExcludeArr{i+1};
    tempSub = subsystemsToExcludeArr{i+1};
    %length(tempSub('Exchange/demand reaction'))
    excludeArr{i+1} = 1:realI;


    reducedModel2 = makeReducedModel(model,subsystemsToAdd2,subsystemsToExcludeIdxs3);
    reducedModelArr2{i+2} = reducedModel2;
end
if 1
parfor i=0:length(reactionsToAdd)
    realI = length(reactionsToAdd)-i;
    reducedModel2 = reducedModelArr2{i+1};
    
    [falconfluxes falconobj] = runFluxMethod(data_257,ids_257,'test',reducedModel2,'FALCON',stds_257,'BIOMASS');
    ratio = falconfluxes(strcmp(reducedModel2.rxns,'TR_glc_D[c]'))/ ...
         -falconfluxes(strcmp(reducedModel2.rxns,'TR_lac_L[c]'));
    if any(strcmp(reducedModel2.rxns,'EX_glc(e)')) && any(strcmp(reducedModel2.rxns,'EX_lac_L(e)'))
        ratio = ratio - falconfluxes(strcmp(reducedModel2.rxns,'EX_glc(e)'))/ ...
         falconfluxes(strcmp(reducedModel2.rxns,'EX_lac_L(e)'));
    end
    ratioArr2(i+2) = ratio;
    disp(reactionsToAdd{i+1})
    disp(i)
end
end
end
end


    

