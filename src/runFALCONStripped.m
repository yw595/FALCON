function [v_falcon objValue cost_rev corr_rho rxn_exp_md_rev] = runFALCONStripped(model,expressionIDs,expData,expressionSDs,nReps,overwriteSDs,overwriteMeans)

    if ~exist('overwriteSDs','var')
        overwriteSDs = containers.Map;
    end
    if ~exist('overwriteMeans','var')
        overwriteMeans = containers.Map;
    end
    
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
    genedata_filename = '/mnt/vdb/home/ubuntu2/temp.csv';
    writeData({expressionIDs,expData,expressionSDs},genedata_filename,'\t');
    [rxn_exp_md, rxn_exp_sd_md, rxn_rule_group] = ... 
        computeMinDisj(modelIrrev, genedata_filename,0,true,overwriteSDs,overwriteMeans);
    minDisjTime = toc(expTime);
    minDisjTime = 0;

    [v_falconIrr, ~, ~, ~, ~, ~, v_falconIrr_s, ~,~,~,~,fOpt, f_easyLP, v_easyLP, cost_irrev] = falconMulti(modelIrrev, nReps, rxn_exp_md, rxn_exp_sd_md, rxn_rule_group, 'rc', regC, 'minFit', minFit, 'EXPCON', expCon,'FDEBUG',FDEBUG);
    v_falcon = convertIrrevFluxDistribution(v_falconIrr, matchRev);
    v_falcon_s = convertIrrevFluxDistribution(v_falconIrr_s, matchRev);
    validIdxs = ~isnan(rxn_exp_md);
    objValue = fOpt;
    objValueOld = sum(abs(v_falconIrr(validIdxs)-rxn_exp_md(validIdxs))./rxn_exp_sd_md(validIdxs));
    fluxSolSum = sum(abs(v_falconIrr(validIdxs)));
    rxn_exp_md_corr = rxn_exp_md(~isnan(rxn_exp_md));
    v_falconIrr_corr = abs(v_falconIrr(~isnan(rxn_exp_md)));
    if length(rxn_exp_md_corr)>0
        corr_rho = corr(rxn_exp_md_corr,v_falconIrr_corr,'type', 'Spearman');
    else
        corr_rho = 0;
    end
    
    rxn_exp_md_rev = zeros(max(irrev2rev),1);
    for i=1:length(rxn_exp_md_rev)
        rxn_exp_md_rev(i) = max(rxn_exp_md(rev2irrev{i}));
    end
    %for i=1:length(irrev2rev)
    %    rxn_exp_md_rev(irrev2rev(i)) = rxn_exp_md_rev(irrev2rev(i)) + rxn_exp_md(i);
    %end
    %cost_irrev = columnVector(f_easyLP).*columnVector(v_easyLP);
    cost_rev = zeros(max(irrev2rev),1);
    for i=1:length(irrev2rev)
        cost_rev(irrev2rev(i)) = cost_rev(irrev2rev(i)) + cost_irrev(i);
    end
end
