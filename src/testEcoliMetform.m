iAF1260_1 = readCbModel('/home/fs01/yw595/MATLAB/FALCON/input/iAF1260_flux1.xml');

noOxygen = 1;
if noOxygen
    iAF1260_1.lb(strcmp(iAF1260_1.rxns,'EX_o2(e)'))=0;
    iAF1260_1.ub(strcmp(iAF1260_1.rxns,'EX_o2(e)'))=0;
end

iAF1260_1 = changeObjective(iAF1260_1,'Ec_biomass_iAF1260_core_59p81M');
iAF1260_1flux = optimizeCbModel(iAF1260_1);
iAF1260_1flux_biomass = iAF1260_1flux.f;
iAF1260_1flux_dist = iAF1260_1flux.x;

biomassArr = [];
lowerBound = 0;
checkRev = 1;
if checkRev
    lowerBound = -50;
end
for i=30:-5:lowerBound
    ithModel = iAF1260_1;
    for j=1:length(ithModel.rxns)
        if ~isempty(regexp(ithModel.rxnNames{j},'NADH dehydrogenase'))
            if i>=0
                ithModel.ub(j) = i/6;
            else
                ithModel.ub(j) = 0;
                ithModel.lb(j) = i/6;
            end            
        end
    end
    ithModelflux = optimizeCbModel(ithModel);
    ithModelflux_biomass = ithModelflux.f;
    ithModelflux_dist = ithModelflux.x;
    i
    ithModelflux_biomass
    for j=1:length(ithModel.rxns)
        if ~isempty(regexp(ithModel.rxnNames{j},'NADH dehydrogenase'))
            ithModelflux_dist(j)
        end
    end
    biomassArr((i+5-lowerBound)/5) = ithModelflux_biomass;
end





