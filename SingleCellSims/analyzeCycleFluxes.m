function analyzeCycleFluxes(model, G1File, SFile, G2MFile, filterExchAndTrans)
    notExchAndTrans = true(length(model.rxns),1);
    if filterExchAndTrans
    for i=1:length(model.rxns)
    	if(~isempty(regexp( model.subSystems{i},'Transport' )) || ...
	~isempty(regexp( model.subSystems{i},'Exchange' )))
	    notExchAndTrans(i) = false;
	end
    end
    end

    load (G1File); G1Fluxes = v_falcon(notExchAndTrans);
    load (SFile); SFluxes = v_falcon(notExchAndTrans);
    load (G2MFile); G2MFluxes = v_falcon(notExchAndTrans);
    modelRxnNames = cell(length(model.rxns),1);
    for i=1:length(model.rxns)
    	modelRxnNames{i} = [model.rxnNames{i} ' (Abbr: ' model.rxns{i} ')'];
	if(length(modelRxnNames{i})>100)
	    modelRxnNames{i} = ['(Abbr: ' model.rxns{i} ')'];
	end
    end
    modelRxnNames = modelRxnNames(notExchAndTrans);

    G1SDiffs = G1Fluxes - SFluxes;
    SG2MDiffs = SFluxes - G2MFluxes;
    [junk sortIdxs] = sort(abs(G1SDiffs),'descend');
    G1SAbsDiffsG1TopTen = G1Fluxes(sortIdxs(1:10));
    G1SAbsDiffsSTopTen = SFluxes(sortIdxs(1:10));
    G1SAbsDiffsTopTenNames = modelRxnNames(sortIdxs(1:10));
    [junk sortIdxs] = sort(abs(SG2MDiffs),'descend');
    SG2MAbsDiffsSTopTen = SFluxes(sortIdxs(1:10));
    SG2MAbsDiffsG2MTopTen = G2MFluxes(sortIdxs(1:10));
    SG2MAbsDiffsTopTenNames = modelRxnNames(sortIdxs(1:10));

    G1Zero = G1Fluxes == 0; G1Nonzero = ~G1Zero;
    G1FluxesZero = G1Fluxes(G1Zero); G1FluxesNonzero = G1Fluxes(G1Nonzero);
    SFluxesZero1 = SFluxes(G1Zero); SFluxesNonzero1 = SFluxes(G1Nonzero);
    G1SDiffsZero = G1SDiffs(G1Zero); G1SDiffsNonzero = G1SDiffs(G1Nonzero);
    G1SDiffsZeroNames = modelRxnNames(G1Zero); G1SDiffsNonzeroNames = modelRxnNames(G1Nonzero);
    
    SZero = SFluxes == 0; SNonzero = ~SZero;
    SFluxesZero2 = SFluxes(SZero); SFluxesNonzero2 = SFluxes(SNonzero);
    G2MFluxesZero = G2MFluxes(SZero); G2MFluxesNonzero = G2MFluxes(SNonzero);
    SG2MDiffsZero = SG2MDiffs(SZero); SG2MDiffsNonzero = SG2MDiffs(SNonzero);
    SG2MDiffsZeroNames = modelRxnNames(SZero); SG2MDiffsNonzeroNames = modelRxnNames(SNonzero);
    
    [junk sortIdxs] = sort(abs(G1SDiffsZero),'descend');
    G1SRelDiffsZeroG1TopTen = G1FluxesZero(sortIdxs(1:10));
    G1SRelDiffsZeroSTopTen = SFluxesZero1(sortIdxs(1:10));
    G1SRelDiffsZeroTopTenNames = G1SDiffsZeroNames(sortIdxs(1:10));
    [junk sortIdxs] = sort(abs(SG2MDiffsZero),'descend');
    SG2MRelDiffsZeroSTopTen = SFluxesZero2(sortIdxs(1:10));
    SG2MRelDiffsZeroG2MTopTen = G2MFluxesZero(sortIdxs(1:10));
    SG2MRelDiffsZeroTopTenNames = SG2MDiffsZeroNames(sortIdxs(1:10));

    [junk sortIdxs] = sort(abs(G1SDiffsNonzero),'descend');
    G1SRelDiffsNonzeroG1TopTen = G1FluxesNonzero(sortIdxs(1:10));
    G1SRelDiffsNonzeroSTopTen = SFluxesNonzero1(sortIdxs(1:10));
    G1SRelDiffsNonzeroTopTenNames = G1SDiffsNonzeroNames(sortIdxs(1:10));
    [junk sortIdxs] = sort(abs(SG2MDiffsNonzero),'descend');
    SG2MRelDiffsNonzeroSTopTen = SFluxesNonzero2(sortIdxs(1:10));
    SG2MRelDiffsNonzeroG2MTopTen = G2MFluxesNonzero(sortIdxs(1:10));
    SG2MRelDiffsNonzeroTopTenNames = SG2MDiffsNonzeroNames(sortIdxs(1:10));

    save('cycleFluxesExportToPython.mat', ...
    'G1SAbsDiffsG1TopTen','G1SAbsDiffsSTopTen','G1SAbsDiffsTopTenNames', ...
    'SG2MAbsDiffsSTopTen','SG2MAbsDiffsG2MTopTen','SG2MAbsDiffsTopTenNames', ...
    'G1SRelDiffsZeroG1TopTen','G1SRelDiffsZeroSTopTen','G1SRelDiffsZeroTopTenNames',...
    'SG2MRelDiffsZeroSTopTen','SG2MRelDiffsZeroG2MTopTen','SG2MRelDiffsZeroTopTenNames',...
    'G1SRelDiffsNonzeroG1TopTen','G1SRelDiffsNonzeroSTopTen','G1SRelDiffsNonzeroTopTenNames',...
    'SG2MRelDiffsNonzeroSTopTen','SG2MRelDiffsNonzeroG2MTopTen','SG2MRelDiffsNonzeroTopTenNames');
end