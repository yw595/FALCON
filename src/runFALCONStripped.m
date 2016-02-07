function v_falcon = runFALCONStripped(model,expressionIDs,expData,expressionSDs,nReps)
    
    %Whether to use new complexation method in Lee method
    useMinDisj = true;
    expCon = false;
    minFit = 0.0;
    regC = 0;

    nrxns = length(model.rxns);

    [modelIrrev, matchRev, rev2irrev, irrev2rev] = convertToIrreversible(model);

    % map gene weighting to reaction weighting

    expTime = tic;
    genedata_filename = '/home/ubuntu/temp.csv';
    writeData({expressionIDs,expData,expressionSDs},genedata_filename,'\t');
    [rxn_exp_md, rxn_exp_sd_md, rxn_rule_group] = ... 
        computeMinDisj(modelIrrev, genedata_filename);
    minDisjTime = toc(expTime);
    minDisjTime = 0;

    [v_falconIrr, ~, ~, ~, ~, ~, v_falconIrr_s] = falconMulti(modelIrrev, nReps, rxn_exp_md, rxn_exp_sd_md, rxn_rule_group, 'rc', regC, 'minFit', minFit, 'EXPCON', expCon);
    v_falcon = convertIrrevFluxDistribution(v_falconIrr, matchRev);
    v_falcon_s = convertIrrevFluxDistribution(v_falconIrr_s, matchRev);
    disp(sum(v_falcon))
end