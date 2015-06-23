function [absDiffs1TopTen absDiffs2TopTen absDiffsTopTenNames ...
	 relDiffsZero1TopTen relDiffsZero2TopTen relDiffsZeroTopTenNames ...
	 relDiffsNonzero1TopTen relDiffsNonzero2TopTen relDiffsNonzeroTopTenNames ...
	 SSAbsDiffs1TopTen SSAbsDiffs2TopTen SSAbsDiffsTopTenNames ...
	 SSRelDiffsZero1TopTen SSRelDiffsZero2TopTen SSRelDiffsZeroTopTenNames ...
	 SSRelDiffsNonzero1TopTen SSRelDiffsNonzero2TopTen SSRelDiffsNonzeroTopTenNames] ...
	 = analyzeDiffFluxes(model, fluxFile1, fluxFile2, filterExchAndTrans)
    notExchAndTrans = true(length(model.rxns),1);
    if filterExchAndTrans
    for i=1:length(model.rxns)
    	if(~isempty(regexp( model.subSystems{i},'Transport' )) || ...
	~isempty(regexp( model.subSystems{i},'Exchange' )) || ...
	strcmp( model.subSystems{i},'' ))
	    notExchAndTrans(i) = false;
	end
    end
    end

    load (fluxFile1); fluxes1 = v_falcon(notExchAndTrans);
    load (fluxFile2); fluxes2 = v_falcon(notExchAndTrans);
    modelRxnNames = cell(length(model.rxns),1);
    for i=1:length(model.rxns)
    	modelRxnNames{i} = [model.rxnNames{i} ' (Abbr: ' model.rxns{i} ') (Subsystem: ' model.subSystems{i} ')'];
	if(length(modelRxnNames{i})>100)
	    modelRxnNames{i} = ['(Abbr: ' model.rxns{i} ') (Subsystem: ' model.subSystems{i} ')'];
	end
    end
    modelRxnNames = modelRxnNames(notExchAndTrans);
    modelSubsystems = model.subSystems;
    modelSubsystems = modelSubsystems(notExchAndTrans);    

    [absDiffs1TopTen absDiffs2TopTen absDiffsTopTenNames ...
	 relDiffsZero1TopTen relDiffsZero2TopTen relDiffsZeroTopTenNames ...
	 relDiffsNonzero1TopTen relDiffsNonzero2TopTen relDiffsNonzeroTopTenNames] ...
	 = calcDiffs(fluxes1, fluxes2, modelRxnNames);

    uniqueSS = unique(model.subSystems);
    uniqueSSTemp = {};
    for i=1:length(uniqueSS)
    	if(~isempty(regexp( uniqueSS{i},'Transport' )) || ...
	~isempty(regexp( uniqueSS{i},'Exchange' )) || ...
	strcmp( uniqueSS{i},'' ))
	    ;
	else
	    uniqueSSTemp{end+1} = uniqueSS{i};
	end
    end
    uniqueSS = uniqueSSTemp;
    SSfluxes1 = []; SSfluxes2 = [];
    for i=1:length(uniqueSS)
    	SSfluxes1(i) = mean(fluxes1(strcmp( modelSubsystems,uniqueSS{i} )));
	SSfluxes2(i) = mean(fluxes2(strcmp( modelSubsystems,uniqueSS{i} )));
    end

    SSfluxes1 = SSfluxes1'; SSfluxes2 = SSfluxes2'; uniqueSS = uniqueSS';
    [SSAbsDiffs1TopTen SSAbsDiffs2TopTen SSAbsDiffsTopTenNames ...
	 SSRelDiffsZero1TopTen SSRelDiffsZero2TopTen SSRelDiffsZeroTopTenNames ...
	 SSRelDiffsNonzero1TopTen SSRelDiffsNonzero2TopTen SSRelDiffsNonzeroTopTenNames] ...
	 = calcDiffs(SSfluxes1, SSfluxes2, uniqueSS);
end