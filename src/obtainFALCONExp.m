function expr = obtainFALCONExp(expressionData,expressionIDs,model,expressionSDs)

    overwriteSDs = containers.Map;
    overwriteMeans = containers.Map;
    [modelIrrev, matchRev, rev2irrev, irrev2rev] = convertToIrreversible(model);
    genedata_filename = '/mnt/vdb/home/ubuntu2/temp.csv';
    writeData({expressionIDs,expressionData,expressionSDs},genedata_filename,'\t');
    [rxn_exp_md, rxn_exp_sd_md, rxn_rule_group] = ... 
        computeMinDisj(modelIrrev, genedata_filename,0,true,overwriteSDs,overwriteMeans);
    expr = convertIrrevFluxDistribution(rxn_exp_md, matchRev);
