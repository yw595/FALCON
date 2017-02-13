ecoliModel = readCbModel('/home/fs01/yw595/MATLAB/FALCON/input/iAF1260_flux1.xml');
rxnsNum = length(ecoliModel.rxns);
rxnsInter = length(intersect(origRecon2.rxns,ecoliModel.rxns));
[~,~,ib1] = intersect(origRecon2.rxns,ecoliModel.rxns);
[~,~,ib2] = intersect(origRecon2.rxnNames,ecoliModel.rxnNames);
rxnsAllInter = length(union(ib1,ib2));
metsNum = length(ecoliModel.mets);
metsInter = length(intersect(origRecon2.mets,ecoliModel.mets));
[~,~,ib1] = intersect(origRecon2.mets,ecoliModel.mets);
[~,~,ib2] = intersect(origRecon2.metNames,ecoliModel.metNames);
metsAllInter = length(union(ib1,ib2));

xvals=[1,1,1,2,2,2];yvals=[rxnsNum-rxnsAllInter-rxnsInter,rxnsAllInter-rxnsInter,rxnsInter,metsNum-metsAllInter-metsInter,metsAllInter-metsInter,metsInter];
xlabels={};grouplabels={};
for i=1:length(xvals)
    if xvals(i)==1
        xlabels{i} = 'Reactions';
    else
        xlabels{i} = 'Metabolites';
    end
    if mod(i,3)==1
        grouplabels{i} = 'Total Number of Entries';
    elseif mod(i,3)==2
        grouplabels{i} = 'Number of Overlaps Based on Either Long or Short ID Match';
    elseif mod(i,3)==0
        grouplabels{i} = 'Number of Overlaps Based on Short ID Match';
    end
end

writeData({xvals,yvals,xlabels,grouplabels},'/home/fs01/yw595/ecoliMap.txt','\t',{'xvals','yvals','xlabels','grouplabels'});










