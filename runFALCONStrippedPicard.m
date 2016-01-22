function v_falcon = runFALCONStrippedPicard(model, genedata_filename, nReps, outputDir)
    
    %Whether to use new complexation method in Lee method
    useMinDisj = true;
    expCon = false;
    minFit = 0.0;
    regC = 0;

    nrxns = length(model.rxns);

    [modelIrrev, matchRev, rev2irrev, irrev2rev] = convertToIrreversible(model);

    % load transcript data
    %genedata	= importdata(genedata_filename);
    %genenames	= genedata.textdata(:,1);
    %genenames(1)= [];
    %gene_exp	= genedata.data(:,1);
    %gene_exp_sd	= genedata.data(:,2);
    FI = fopen(genedata_filename);
    dataFields = textscan(FI,'%s\t%f\t%f','Delimiter','\t');
    fclose(FI);
    genenames = dataFields{1};
    gene_exp = dataFields{2};
    gene_exp_sd = dataFields{3};

    % map gene weighting to reaction weighting

    expTime = tic;
    [rxn_exp_md, rxn_exp_sd_md, rxn_rule_group] = ... 
        computeMinDisj(modelIrrev, genedata_filename);
    %rxn_exp_md
    minDisjTime = toc(expTime);
    minDisjTime = 0;
    
    if ~exist(outputDir,'dir')
        system(['mkdir ' outputDir]);
    end

    [~, genedata_filename_name, genedata_filename_ext] = fileparts(genedata_filename);
    stripped_genedata_filename = [genedata_filename_name genedata_filename_ext];
    [v_falconIrr, ~, ~, ~, ~, ~, v_falconIrr_s] =         ...
        falconMulti(modelIrrev, nReps, rxn_exp_md,        ...
        rxn_exp_sd_md, rxn_rule_group, 'rc', regC,        ...
        'minFit', minFit, 'EXPCON', expCon, outputDir, 	  ...
	stripped_genedata_filename);
    v_falcon = convertIrrevFluxDistribution(v_falconIrr, matchRev);
    v_falcon_s = convertIrrevFluxDistribution(v_falconIrr_s, matchRev);
    disp(sum(v_falcon))
    
    save([outputDir filesep stripped_genedata_filename '_falcon_flux.mat'], 'v_falcon');
    
    fluxOutputFI = fopen([outputDir filesep stripped_genedata_filename '.flux'],'w');
    for i=1:length(v_falcon)
    	fprintf(fluxOutputFI,['R_' model.rxns{i} ',' num2str(v_falcon(i)) '\n']);
    end
    fclose(fluxOutputFI);
end