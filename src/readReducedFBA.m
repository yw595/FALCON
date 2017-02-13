transferFiles = dir('/home/fs01/yw595/transferTo');
%transferFiles = transferFiles.name;
count = 0;
subsystemsArr = {};
cutoffsArr = [];
subCutoffsArr = [];
uniqueSubs = strrep(strrep(unique(origRecon2.subSystems),' ','_'),',','_');
for i=1:length(transferFiles)
    %disp(transferFiles(i).name);
    if count>10
        %break
    end
    prefix = 'reducedFBASubAllDiff';
    if ~isempty(regexp(transferFiles(i).name,'reducedFBASubAll')) && isempty(regexp(transferFiles(i).name,'Fine')) && ~isempty(regexp(transferFiles(i).name,'Diff'))
    count=count+1;
    %disp(transferFiles(i).name)
    subsystem = transferFiles(i).name;
    subsystem = subsystem(length(prefix)+1:end-4);
    while ~isempty(str2num(subsystem(end))) && length(subsystem)>=2
        subsystem = subsystem(1:end-1);
    end
    if ~any(strcmp(uniqueSubs,subsystem))
    FI = fopen(['/home/fs01/yw595/transferTo' filesep transferFiles(i).name]);
    dataFields = textscan(FI,'%s%d%s%s%s%f%s%s','Delimiter','\t','Headerlines',0);
    transferRxns = dataFields{1};
    transferFlux = abs(dataFields{6});
    [uniqueRxns,~,uniqueIdxs] = unique(transferRxns);
    uniqueFlux = [];
    for j=1:max(uniqueIdxs)
        uniqueFlux(j) = sum(transferFlux(uniqueIdxs==j));
    end
    [uniqueFlux,sortIdxs] = sort(uniqueFlux,'descend');
    uniqueRxns = uniqueRxns(sortIdxs);
    cutoffIdx = min(find(cumsum(uniqueFlux) >= .9*sum(uniqueFlux)));
    cutRxns = uniqueRxns(1:cutoffIdx);
    uniqueSubsystems = {};%cellfun(@(x) origRecon2.subSystems{strcmp(origRecon2.rxns,x)}, cutRxns,'UniformOutput',0);
    uniqueSubsystemsCounts = containers.Map;
    for j=1:length(cutRxns)
        isSubsystem = find(strcmp(origRecon2.rxns,cutRxns{j}));
        if ~isempty(isSubsystem)
            uniqueSubsystem = origRecon2.subSystems{strcmp(origRecon2.rxns,cutRxns{j})};
            uniqueSubsystems{j} = uniqueSubsystem;
        else
            uniqueSubsystems{j} = '';
        end

        if ~isKey(uniqueSubsystemsCounts,uniqueSubsystems{j})
            uniqueSubsystemsCounts(uniqueSubsystems{j}) = 1;
        else
            uniqueSubsystemsCounts(uniqueSubsystems{j}) = uniqueSubsystemsCounts(uniqueSubsystems{j}) + 1;
        end
    end
    %nonsense = nonsense+1;
    disp(subsystem)
    disp(count)
    disp(cutoffIdx)
    cutoffsArr(count) = cutoffIdx;
    subCutoffsArr(count) = length(keys(uniqueSubsystemsCounts));
    subsystemsArr{count} = subsystem;
    fclose(FI);
    %nonsense = nonsense+1;
    end
    end
end







    



                