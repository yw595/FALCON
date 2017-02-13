[numarray textarray raw] = xlsread('/home/fs01/yw595/NMT113-gene_list.xlsx');
[intersectGenes, ~, intersectIdxs] = intersect(textarray(:,1),yeastModel.genes);
intersectRxnNames = yeastModel.rxnNames(sum(yeastModel.rxnGeneMat(:,intersectIdxs),2)~=0);
intersectRxns = yeastModel.rxns(sum(yeastModel.rxnGeneMat(:,intersectIdxs),2)~=0);
multimapGenes = cellfun(@(x) yeastModel.genes(find(yeastModel.rxnGeneMat(strcmp(yeastModel.rxns,x),:))),intersectRxns,'UniformOutput',0);
multimapRxns = cellfun(@(x) yeastModel.rxns(find(yeastModel.rxnGeneMat(:,strcmp(yeastModel.genes,x)))),intersectGenes,'UniformOutput',0);

multimapRxnNames = {};
for i=1:length(multimapRxns)
    multimapRxnNames{i} = {};
    for j=1:length(multimapRxns{i})
        multimapRxnNames{i}{j} = yeastModel.rxnNames{strcmp(yeastModel.rxns,multimapRxns{i}{j})};
    end
end
for i=1:length(multimapRxnNames)
    disp(intersectGenes{i})
    dispString = multimapRxnNames{i}{1};
    for j=2:length(multimapRxnNames{i})
        dispString = [dispString ' : ' multimapRxnNames{i}{j}];
    end
    disp(dispString)
end

multimapGeneDesc = {};
for i=1:length(multimapGenes)
    multimapGeneDesc{i} = {};
    for j=1:length(multimapGenes{i})
        if any(strcmp(intersectGenes,multimapGenes{i}{j}))
            multimapGeneDesc{i}{j} = textarray{strcmp(textarray(:,1),multimapGenes{i}{j}),3};
        else
            multimapGeneDesc{i}{j} = 'None';
        end
    end
end
for i=1:length(multimapGeneDesc)
    disp(intersectRxns{i})
    dispString = multimapGeneDesc{i}{1};
    for j=2:length(multimapGeneDesc{i})
        dispString = [dispString ' : ' multimapGeneDesc{i}{j}];
    end
    disp(dispString)
end













