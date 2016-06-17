function printModel(model,fluxDist,outputFile,frequentFile)

FI = fopen(outputFile,'w');
for i=1:length(model.rxns)
    metIdxs = find(model.S(:,i));
    for j=1:length(metIdxs)
        fprintf(FI,'%s\t%d\t%s\t%s\t%s\t%f\t%s\t%s',model.rxns{i},full(model.S(metIdxs(j),i)),model.mets{metIdxs(j)},model.subSystems{i},model.metKeggID{metIdxs(j)},fluxDist(i),model.rxnNames{i},model.metNames{metIdxs(j)});
        if i ~= length(model.rxns) || j~=length(metIdxs)
            fprintf(FI,'\n');
        end
    end
end
fclose(FI);

frequentFI = fopen(frequentFile,'w');
for i=1:length(model.mets)
    numRxnsConn = sum(abs(model.S(i,:)));
    if numRxnsConn > 10
        fprintf(frequentFI,'%s\n',model.mets{i});
    end
end
fclose(frequentFI);