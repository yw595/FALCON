function [ratioArr distArr corrArr avgCostArr costArr expArr] = runFALCONModelArr(modelArr,skipFirst,dispIters,option)

inputDirmRNA = '/home/fs01/yw595/MATLAB/FALCON/input/NCI60Sims/nci60mRNA';
[ids_257 data_257 stds_257] = readExpressionFile([inputDirmRNA '/UACC_257.csv']);

ratioArr = [];
distArr = {};
corrArr = [];
avgCostArr = [];
costArr = {};
expArr = {};
startIdx = 1;
if skipFirst
    startIdx = 2;
end
parfor i=1:length(modelArr)
    overwriteSDs = containers.Map;
    overwriteMeans = containers.Map;
    %if exist('option','var')
        if strcmp(option,'cents')
            cents = betweenness_centrality(sparse(makeRxnConnMatrix(modelArr{i},1)));
            cents1 = cents;
            cents1(cents1==0) = min(cents1(cents1~=0))/10;
            for j=1:length(modelArr{i}.rxns)
                overwriteSDs(num2str(j)) = 1/cents1(j);
            end
        elseif strcmp(option,'zero')
            [modelIrrev, matchRev, rev2irrev, irrev2rev] = convertToIrreversible(modelArr{i});
            genedata_filename = '/home/fs01/yw595/temp.csv';
            [rxn_exp_md, rxn_exp_sd_md, rxn_rule_group] = computeMinDisj(modelIrrev, genedata_filename);
            rxn_exp_md_rev = zeros(length(modelArr{i}.rxns),1);
            for j=1:length(irrev2rev)
                rxn_exp_md_rev(irrev2rev(j)) = rxn_exp_md_rev(irrev2rev(j)) + rxn_exp_md(j);
            end
            for j=1:length(modelArr{i}.rxns)
                if rxn_exp_md_rev(j)==0
                    overwriteSDs(num2str(j)) = .001;
                else
                    overwriteSDs(num2str(j)) = 1;
                end
            end
        end
        %end
    
    initCobraToolbox;
    [falconfluxes falconobj cost_rev corr_rho rxn_exp_md_rev] = runFluxMethod(data_257,ids_257,'test',modelArr{i},'FALCON',stds_257,'BIOMASS',overwriteSDs,overwriteMeans);
    ratioArr(i) = falconfluxes(strcmp(modelArr{i}.rxns,'TR_glc_D[c]'))/ ...
         -falconfluxes(strcmp(modelArr{i}.rxns,'TR_lac_L[c]'));
    if any(strcmp(modelArr{i}.rxns,'EX_glc(e)')) && any(strcmp(modelArr{i}.rxns,'EX_lac_L(e)'))
        ratioArr(i) = ratioArr(i) - falconfluxes(strcmp(modelArr{i}.rxns,'EX_glc(e)'))/ ...
         falconfluxes(strcmp(modelArr{i}.rxns,'EX_lac_L(e)'));
    end
    distArr{i} = falconfluxes;
    corrArr(i) = corr_rho;
    avgCostArr(i) = sum(cost_rev)/length(cost_rev);
    costArr{i} = cost_rev;
    expArr{i} = rxn_exp_md_rev;
    if dispIters
        disp(i);
    end
end

end