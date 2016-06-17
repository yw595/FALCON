function [v_falcon objValue] = runFALCONStripped(model,expressionIDs,expData,expressionSDs,nReps)
    
    %Whether to use new complexation method in Lee method
    useMinDisj = true;
    expCon = false;
    minFit = 0.0;
    regC = 0;
    FDEBUG = false;

    nrxns = length(model.rxns);

    [modelIrrev, matchRev, rev2irrev, irrev2rev] = convertToIrreversible(model);

    % map gene weighting to reaction weighting

    expTime = tic;
    genedata_filename = '/home/fs01/yw595/temp.csv';
    writeData({expressionIDs,expData,expressionSDs},genedata_filename,'\t');
    [rxn_exp_md, rxn_exp_sd_md, rxn_rule_group] = ... 
        computeMinDisj(modelIrrev, genedata_filename);
    minDisjTime = toc(expTime);
    minDisjTime = 0;

    [v_falconIrr, ~, ~, ~, ~, ~, v_falconIrr_s, ~,~,~,~,fOpt] = falconMulti(modelIrrev, nReps, rxn_exp_md, rxn_exp_sd_md, rxn_rule_group, 'rc', regC, 'minFit', minFit, 'EXPCON', expCon,'FDEBUG',FDEBUG);
    disp(v_falconIrr)
    v_falcon = convertIrrevFluxDistribution(v_falconIrr, matchRev);
    disp(v_falcon)
    v_falcon_s = convertIrrevFluxDistribution(v_falconIrr_s, matchRev);
    validIdxs = ~isnan(rxn_exp_md);
    objValue = fOpt;
    objValueOld = sum(abs(v_falconIrr(validIdxs)-rxn_exp_md(validIdxs))./rxn_exp_sd_md(validIdxs));
    fluxSolSum = sum(abs(v_falconIrr(validIdxs)));
    disp(objValue)
    disp(objValueOld)
    disp(fluxSolSum)
end