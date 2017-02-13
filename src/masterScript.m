model = origRecon2;
paramSets = {};
inputDirmRNA = '/home/fs01/yw595/MATLAB/FALCON/input/NCI60Sims/nci60mRNA';
transferDir = '/home/fs01/yw595/MATLAB/FALCON/transfer';
[ids_257 data_257 stds_257] = readExpressionFile([inputDirmRNA '/UACC_257.csv']);

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

if 0 %make nanZeroExp data
[modelIrrev, matchRev, rev2irrev, irrev2rev] = convertToIrreversible(origRecon2);
[rxn_exp_md, rxn_exp_sd_md, rxn_rule_group] = computeMinDisj(modelIrrev, genedata_filename);
rxn_exp_md_rev = zeros(length(origRecon2.rxns),1);
for i=1:length(irrev2rev)
    rxn_exp_md_rev(irrev2rev(i)) = rxn_exp_md_rev(irrev2rev(i)) + rxn_exp_md(i);
end
subNanExpPcts = []; subZeroExpPcts = []; uniqSubs = unique(origRecon2.subSystems);
for i=1:length(uniqSubs)
    subExp = rxn_exp_md_rev(strcmp(origRecon2.subSystems,uniqSubs{i}));
    subNanExpPcts(i) = sum(isnan(subExp))/length(subExp);
    subZeroExpPcts(i) = sum(subExp==0)/length(subExp);
end
uniqSubsShow = cell(length(uniqSubs)*2,1);
uniqSubsShow(1:2:length(uniqSubsShow)-1) = cellfun(@(x) [x ' nan'], uniqSubs,'UniformOutput',0);
uniqSubsShow(2:2:length(uniqSubsShow)) = cellfun(@(x) [x ' zero'], uniqSubs,'UniformOutput',0);
subPctsShow = zeros(length(uniqSubs)*2,1);
subPctsShow(1:2:length(subPctsShow)-1) = subNanExpPcts;
subPctsShow(2:2:length(subPctsShow)) = subZeroExpPcts;
groupFill = cell(length(uniqSubs)*2,1);
for i=1:length(groupFill)
    if mod(i,2)==1
        groupFill{i} = 'nan';
    else
        groupFill{i} = 'zero';
    end
end
writeData({uniqSubsShow,subPctsShow,groupFill},'/home/fs01/yw595/nanZeroExpPcts.txt','\t',{'subsystems','pcts','groupFill'});
end



if 0 %make initial curated reducedModel
subsystemsToExcludeIdxs = containers.Map;
subsystemsToAdd = {};
reducedModel = makeReducedModel(model,subsystemsToAdd,subsystemsToExcludeIdxs);
end

if 0 %make reducedModelArr and subsystemsAddedArr
reducedModelArr = {reducedModel};
subsystemsAddedArr = {'None'};
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

while length(unique(reducedModelArr{end}.subSystems)) ~= length(unique(origRecon2.subSystems))
    mostComm = getMostCommSub(reducedModelArr{end},origRecon2);
    mostComm
    
    subsystemsToAdd{end+1} = mostComm;
    subsystemsToExcludeIdxs(mostComm) = [];
    reducedModelArr{end+1} = makeReducedModel(model,subsystemsToAdd,subsystemsToExcludeIdxs);
    subsystemsAddedArr{end+1} = mostComm;    
end

for i=1:length(subsystemsAddedArr)
    subsystemsAddedArr{i} = strrep(subsystemsAddedArr{i},'/','_');
    subsystemsAddedArr{i} = strrep(subsystemsAddedArr{i},' ','_');
    subsystemsAddedArr{i} = strrep(subsystemsAddedArr{i},',','_');
end
end

%[falconfluxes falconobj cost_rev corr_rho rxn_exp_md_rev] = runFluxMethod(data_257,ids_257,'test',reducedModel,'FALCON',stds_257,'BIOMASS');

%writeData({subsystemsAddedArr(2:end)},[transferDir '/subsystemsAddedArr.txt'],'\t');

if 0 %run FALCON for each subsystem model
[ratioArr distArr corrArr avgCostArr costArr expArr] = runFALCONModelArr(reducedModelArr,0,0)
end

if 0 %write subsystemCorrAndCost data
    writeData({subsystemsAddedArr,corrArr,avgCostArr},[transferDir '/subsystemCorrAndCost.txt'],'\t',{'subsystem','corr_with_exp','avg_cost'});
end

if 0 %write out scatter text files between expression and flux for
     %subsystem added models
for i=1:length(reducedModelArr)
    writeOut1Arr = expArr{i}; writeOut2Arr = distArr{i};
    nanidxs = isnan(writeOut1Arr) | isnan(writeOut2Arr);
    writeOut1Arr = writeOut1Arr(~nanidxs); writeOut2Arr = writeOut2Arr(~nanidxs); 
    %writeOut1Arr(mod(writeOut1Arr,1)==0) = writeOut1Arr(mod(writeOut1Arr,1)==0) + eps;
    %writeOut2Arr(mod(writeOut2Arr,1)==0) = writeOut2Arr(mod(writeOut2Arr,1)==0) + eps;
    %writeOut1 = arrayfun(@(x) x(1)+eps, writeOut1Arr, 'UniformOutput', 0);
    %writeOut2 = arrayfun(@(x) abs(x(1))+eps, writeOut2Arr, 'UniformOutput', 0);
    %writeOut1(2:end+1) = writeOut1(1:end);
    %writeOut2(2:end+1) = writeOut2(1:end);
    %writeOut1{1} = 'exp'; writeOut2{1} = 'flux';
    writeData({writeOut1Arr,writeOut2Arr},[transferDir '/scatterSubAll' num2str(i) subsystemsAddedArr{i} '.txt'],'\t',{'exp','flux'});
end
end

if 1 %randomization trials on subsystem models
single = 0;
iterations = 5;
if single
    iterations = 1;
end

reducedConstrainedModelArrTemp = reducedConstrainedModelArr;
ratioArrArr = {};
for z=1:iterations
    if 1 % run FALCON on finegrained reducedModelArr2
        [ratioArr distArr corrArr avgCostArr costArr expArr] = runFALCONModelArr(reducedConstrainedModelArrTemp,0,0,'zero');
    end
    ratioArrArr{z} = ratioArr;
end

ratioArrArrWrite = [];
for z=1:iterations
    for z1=1:length(reducedConstrainedModelArrTemp)
        ratioArrArrWrite(z1,z) = ratioArrArr{z}(z1);
    end
end

ratioArrArrWrite(2,4) = nan;

medianWrite = [];
lowerWrite = [];
upperWrite = [];
nanInfCounts = [];
for z1=1:length(reducedConstrainedModelArrTemp)
    validRatios = ratioArrArrWrite(z1,~(isinf(ratioArrArrWrite(z1,:)) | isnan(ratioArrArrWrite(z1,:))));
    if ~isempty(validRatios)
        medianWrite(z1) = median(validRatios);
        lowerWrite(z1) = min(validRatios);
        upperWrite(z1) = max(validRatios);
    else
        medianWrite(z1) = 0;
        lowerWrite(z1) = 0;
        upperWrite(z1) = 0;
    end
    nanInfCounts(z1) = sum(isinf(ratioArrArrWrite(z1,:)) | isnan(ratioArrArrWrite(z1,:)));
end

writeData({subsystemsAddedArr,medianWrite,lowerWrite,upperWrite,nanInfCounts},'/home/fs01/yw595/GlcLacRatioZeroIter.txt','\t',{'subsystem','flux','lower','upper','naninf'});
end

if 0 %randomization trials on fine-grained models
single = 1;
randomize = 0;
iterations = 5;
if single
    iterations = 1;
end

reactionsAddedArr2Arr = {};
ratioArr2Arr = {};
for z=1:iterations
    [reducedModelArr2,reactionsAddedArr2,ratioArr2,distArr2,corrArr2,avgCostArr2,costArr2,expArr2,subsystemsToAdd2] = makeFineGrainedReducedModelArr(76,randomize,model,reducedModelArr,ratioArr,distArr,corrArr,avgCostArr,costArr,expArr,subsystemsToAdd);

    if 0 % run FALCON on finegrained reducedModelArr2
        [ratioArr2 distArr2 corrArr2 avgCostArr2 costArr2 expArr2] = runFALCONModelArr(reducedModelArr2,1,1,'default')
    end

    reactionsAddedArr2Arr{z} = reactionsAddedArr2;
    ratioArr2Arr{z} = ratioArr2;
end
end

if 0 %print out flux distributions for subsystem added models
parfor i=1:length(reducedModelArr)
    disp(['Subsystem flux ' num2str(i)])
    printModel(reducedModelArr{i},distArr{i}, ...
               [transferDir '/reducedFBASubAll' subsystemsAddedArr{i} num2str(i) '.txt'], ...
               [transferDir '/frequentMetsReducedSubAll' subsystemsAddedArr{i} num2str(i) '.txt']);
end
end

if 0
parfor i=1:length(reducedModelArr2)
    disp(['Exchange flux ' num2str(i)])
    printModel(reducedModelArr2{i},distArr2{i}, ...
               [transferDir '/reducedFBASubAllFine' reactionsAddedArr2{i} num2str(i) '.txt'], ...
               [transferDir '/frequentMetsReducedSubAllFine' reactionsAddedArr2{i} num2str(i) '.txt']);
end
end

if 0 %calculate and print out differential flux between subsystems
diffDistConstrainedArr = {};
for i=2:length(distConstrainedArr)
    commLength = min(length(distConstrainedArr{i}),length(distConstrainedArr{i-1}));
    diffDistConstrainedArr{i-1} = distConstrainedArr{i}(1:commLength) - distConstrainedArr{i-1}(1:commLength);
end
parfor i=2:length(reducedConstrainedModelArr)
    %disp(['Subsystem diff flux ' num2str(i)])
    %printModel(reducedModelArr{i-1},diffDistArr{i-1}, ...
    %[transferDir '/reducedFBASubAllDiff' subsystemsAddedArr{i} num2str(i) '.txt'], ...
    %[transferDir '/frequentMetsReducedSubAllDiff' subsystemsAddedArr{i} num2str(i) '.txt']);
end
end

if 0 %calculate and print out differential flux between exchange reactions
diffDistArr2 = {};
for i=2:length(distArr2)
    commLength = min(length(distArr2{i}),length(distArr2{i-1}));
    diffDistArr2{i-1} = distArr2{i}(1:commLength) - distArr2{i-1}(1:commLength);
end
parfor i=2:length(reducedModelArr2)
    disp(['Exchange diff flux ' num2str(i)])
    printModel(reducedModelArr2{i-1},diffDistArr2{i-1}, ...
               [transferDir '/reducedFBASubAllDiff' reactionsAddedArr2{i} num2str(i) '.txt'], ...
               [transferDir '/frequentMetsReducedSubAllDiff' reactionsAddedArr2{i} num2str(i) '.txt']);
end
end

if 0
nanInfBins = histc(find(isnan(ratioArr2) | isinf(ratioArr2)),[0,100,200,300,400,500,600,700,Inf]);
nanInfLabels = {}; nanInfPercents = [];
for i=1:length(nanInfBins)-2
    nanInfLabels{i} = [num2str(100*(i-1)+1) '-' num2str(100*i)];
    nanInfPercents(i) = nanInfBins(i)/100;
end
nanInfLabels{end+1} = [num2str(100*(length(nanInfLabels))+1) '-' num2str(length(ratioArr2))];
nanInfPercents(end+1) = nanInfBins(end-1)/mod(length(ratioArr2),100);
writeData({nanInfLabels,nanInfPercents},[transferDir '/nanInfBins.txt'],'\t',{'reactionNumber','percentage'});

writeData({reactionsAddedArr2,ratioArr2},[transferDir '/GlcLacRatioFine.txt'],'\t',{'reactionAdded','ratio'});
end

if 0 % make constrained version of reducedModelArr
reducedConstrainedModelArr = {};
for i=1:length(reducedModelArr)
    reducedConstrainedModelArr{i} = constrainMediumExc(reducedModelArr{i});
end
end

if 0 % make constrained CORE version of reducedModelArr
reducedConstrainedCOREModelArr = {};
[cellLinesArray metsArray coreTable FVAVminArray FVAVmaxArray] = readJainTable();
for i=1:length(reducedModelArr)
    for j=1:length(metsArray)
        disp([num2str(i) ' ' num2str(j)])
        reducedConstrainedCOREModelArr{i,j} = constrainMediumExc(reducedModelArr{i},1,'constrainCORE',{metsArray{j}});
    end
end
end

if 0 % run FALCON with centrality weighting
[ratioConstrainedCentsArr ratioConstrainedCentsArr corrConstrainedCentsArr avgCostConstrainedCentsArr costConstrainedCentsArr expConstrainedCentsArr] = runFALCONModelArr(reducedConstrainedModelArr,0,0,'cents')

writeData({subsystemsAddedArr,ratioConstrainedCentsArr},'/home/fs01/yw595/GlcLacRatioConstrainedCents.txt','\t',{'subsystemAdded','ratio'});
end

if 0 % run FALCON with zero-expression weighting
[ratioConstrainedZeroArr ratioConstrainedZeroArr corrConstrainedZeroArr avgCostConstrainedZeroArr costConstrainedZeroArr expConstrainedZeroArr] = runFALCONModelArr(reducedConstrainedModelArr,0,0,'zero')

writeData({subsystemsAddedArr,ratioConstrainedZeroArr},'/home/fs01/yw595/GlcLacRatioConstrainedZero.txt','\t',{'subsystemAdded','ratio'});
end

if 0 %run FALCON for constrainedCORE models
distConstrainedCOREArr = {};
[cellLinesArray metsArray coreTable FVAVminArray FVAVmaxArray] = readJainTable();
for i=1:length(reducedConstrainedCOREModelArr)
    parfor j=1:length(metsArray)
        initCobraToolbox;
        [falconfluxes falconobj cost_rev corr_rho rxn_exp_md_rev] = runFluxMethod(data_257,ids_257,'test',reducedConstrainedCOREModelArr{i,j},'FALCON',stds_257,'BIOMASS');
        distConstrainedCOREArr{i,j} = falconfluxes;
    end
end
end

if 0 % run FALCON with CORE Medium constraints
[ratioConstrainedArr ratioConstrainedArr corrConstrainedArr avgCostConstrainedArr costConstrainedArr expConstrainedArr] = runFALCONModelArr(reducedConstrainedModelArr,0,0)

writeData({subsystemsAddedArr,ratioConstrainedArr},'/home/fs01/yw595/GlcLacRatioConstrained.txt','\t',{'subsystemAdded','ratio'});
end

useConstrainedCORE = 0;
useConstrainedCents = 0;
useConstrainedZero = 1;
if 0 % calculate exchange flux accuracy
initCobraToolbox;
[cellLinesArray metsArray coreTable FVAVminArray FVAVmaxArray] = readJainTable();
coreCol = coreTable(:,strcmp(cellLinesArray,'UACC_257'));
sensArr = [];
excStatArr = {};
for i=1:length(reducedConstrainedModelArr)
    i
    jainMetsToExcIdxs = loadJainMetsToExcIdxs(metsArray,reducedConstrainedModelArr{i},1);
    FVAVminArray = zeros(length(metsArray),1);
    FVAVmaxArray = zeros(length(metsArray),1);
    for j=1:length(metsArray)
        excIdxs = jainMetsToExcIdxs(metsArray{j});
        if ~isempty(excIdxs)
            j
        end
        for k=1:length(excIdxs)
            if useConstrainedCORE
                FVAModel = changeObjective(reducedConstrainedCOREModelArr{i,j},reducedConstrainedCOREModelArr{i,j}.rxns{excIdxs(k)});
            else
                FVAModel = changeObjective(reducedConstrainedModelArr{i},reducedConstrainedModelArr{i}.rxns{excIdxs(k)});
            end
            FVARes = optimizeCbModel(FVAModel,'min');
            if FVARes.f < FVAVminArray(j)
                FVAVminArray(j) = FVARes.f;
            end
            FVARes = optimizeCbModel(FVAModel,'max');
            if FVARes.f > FVAVmaxArray(j)
                FVAVmaxArray(j) = FVARes.f;
            end
        end
    end

    if useConstrainedCORE
        v_Exc = extractExcFlux(reducedConstrainedCOREModelArr{i,j}, distConstrainedCOREArr{i,j},1);
    elseif useConstrainedCents
        v_Exc = extractExcFlux(reducedConstrainedModelArr{i}, distConstrainedCentsArr{i},1);
    elseif useConstrainedZero
        v_Exc = extractExcFlux(reducedConstrainedModelArr{i}, distConstrainedZeroArr{i},1);
    else
        v_Exc = extractExcFlux(reducedConstrainedModelArr{i}, distConstrainedArr{i},1);
    end
        
    uptakeTruePos=0;
    uptakeFalseNeg=0;
    releaseTruePos=0;
    releaseFalseNeg=0;
    for j=1:length(metsArray)
        disp(metsArray{j})
        if coreCol(j) > 0
            if FVAVmaxArray(j) == 0
                excStatArr{i,j} = 'releaseFVAErr';
            else
                if (v_Exc(j) > 0)
                    disp('releaseTrue')
                    releaseTruePos = releaseTruePos + 1;
                    excStatArr{i,j} = 'releaseTruePos';
                else
                    releaseFalseNeg = releaseFalseNeg + 1;
                    excStatArr{i,j} = 'releaseFalseNeg';
                end
            end
        elseif coreCol(j) < 0
            if FVAVminArray(j) == 0
                excStatArr{i,j} = 'uptakeFVAErr';
            else
                if (v_Exc(j) < 0)
                    disp('uptakeTrue')
                    uptakeTruePos = uptakeTruePos + 1;
                    excStatArr{i,j} = 'uptakeTruePos';
                else
                    uptakeFalseNeg = uptakeFalseNeg + 1;
                    excStatArr{i,j} = 'uptakeFalseNeg';
                end
            end
        else
            excStatArr{i,j} = 'COREZero';
        end
    end
    uptakeFalsePos = releaseFalseNeg;
    uptakeTrueNeg = releaseTruePos;
    releaseFalsePos = uptakeFalseNeg;
    releaseTrueNeg = uptakeTruePos;

    sensArr(i) = (uptakeTruePos + releaseTruePos) / ...
            (uptakeTruePos + uptakeFalseNeg + ...
            releaseTruePos + releaseFalseNeg);
end

end

if 0 % make diff flux heat map data
subsystemsAddedAppendArr = appendNumFunc(subsystemsAddedArr);
x1d = []; y1d = []; diffFluxArr1d = []; xLabel1d = {}; yLabel1d = {};
for i=2:length(reducedConstrainedModelArr)
    subConstrained = clustSubFlux(reducedConstrainedModelArr{i-1},diffDistConstrainedArr{i-1});
    for j=1:length(subsystemsAddedAppendArr)
        x1d(end+1) = j;
        y1d(end+1) = i;
        if isKey(subConstrained,subsystemsAddedArr{j})
            diffFluxArr1d(end+1) = log(subConstrained(subsystemsAddedArr{j}));
        else
            diffFluxArr1d(end+1) = 0;
        end
        xLabel1d{end+1} = subsystemsAddedAppendArr{j};
        yLabel1d{end+1} = subsystemsAddedAppendArr{i};
    end
end

writeData({x1d,y1d,diffFluxArr1d,xLabel1d,yLabel1d},'/home/fs01/yw595/diffFluxHeatmap.txt','\t',{'x','y','diffFlux','sub','model'});
end

if 0 % make exchange statistic heat map data
x1d = []; y1d = []; excStatArr1d = {}; xLabel1d = {}; yLabel1d = {};
for i=1:length(reducedConstrainedModelArr)
    for j=1:length(metsArray)
        x1d(end+1) = j;
        y1d(end+1) = i;
        excStatArr1d{end+1} = excStatArr{i,j};
        xLabel1d{end+1} = metsArray{j};
        yLabel1d{end+1} = subsystemsAddedArr{i};
        appendNum = num2str(i);
        if length(appendNum)<2
            appendNum = ['0' appendNum];
        end
        yLabel1d{end} = [appendNum '_' yLabel1d{end}];
    end
end

if useConstrainedCORE
    writeData({x1d,y1d,excStatArr1d,xLabel1d,yLabel1d},'/home/fs01/yw595/excStatHeatmapConstrainedCORE.txt','\t',{'x','y','excStat','met','model'});
else
    if useConstrainedCents
        writeData({x1d,y1d,excStatArr1d,xLabel1d,yLabel1d},'/home/fs01/yw595/excStatHeatmapConstrainedCents.txt','\t',{'x','y','excStat','met','model'});
    elseif useConstrainedZero
        writeData({x1d,y1d,excStatArr1d,xLabel1d,yLabel1d},'/home/fs01/yw595/excStatHeatmapConstrainedZero.txt','\t',{'x','y','excStat','met','model'});
    else
        writeData({x1d,y1d,excStatArr1d,xLabel1d,yLabel1d},'/home/fs01/yw595/excStatHeatmap.txt','\t',{'x','y','excStat','met','model'});
    end
end
end

writeData({subsystemsAddedArr,ratioConstrainedArrFake},[transferDir '/GlcLacRatioCancer.txt'],'\t',{'subsystemAdded','ratio'});