function minModel = makeMinimalFBAModel(origModel)

minModel = origModel;
keepgoing = 1;
prevZeroIdxs = [];
count = 0;
while keepgoing==1 & count < 500
    count = count+1;
    disp(sum(minModel.lb~=0))
    fbasol = optimizeCbModel(minModel);
    fbafluxes = fbasol.x;
    disp(sum(fbafluxes==0))
    zeroIdxs = find(fbafluxes==0);
    if length(zeroIdxs) > length(prevZeroIdxs)
        prevZeroIdxs = zeroIdxs;
        minModel.lb(fbafluxes==0) = 0;
        minModel.ub(fbafluxes==0) = 0;
    else
        smallestNonzeroIdx = 1;
        absfbafluxes = abs(fbafluxes);
        smallestNonzero = max(absfbafluxes);
        for i=1:length(absfbafluxes)
            if absfbafluxes(i) ~= 0 && absfbafluxes(i) < smallestNonzero
                smallestNonzero = absfbafluxes(i);
                smallestNonzeroIdx = i;
            end
        end
        minModel.lb(smallestNonzeroIdx) = 0;
        minModel.ub(smallestNonzeroIdx) = 0;
        %else
        %keepgoing = 0;
    end
end
%[~ smallestIdx] = min(abs(fbafluxes));
%testMin = minModel; testMin.lb(smallestIdx) = 0; testMin.ub(smallestIdx) = 0;




