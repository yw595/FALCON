function v_Exc = extractExcFlux(model, v_falcon, useReducedNames)

    if ~exist('useReducedNames','var')
        useReducedNames = 0;
    end
    [~, jainMetsArray, ~] = readJainTable();
    jainMetsToExcIdxs = loadJainMetsToExcIdxs(jainMetsArray, model,useReducedNames);
    v_Exc = [];

    for j = 1:length(jainMetsArray)
        jthExcIdxs = jainMetsToExcIdxs(jainMetsArray{j});
	if ~isempty(jthExcIdxs)
            v_Exc(j) = sum(v_falcon(jthExcIdxs));
	else
	    v_Exc(j) = 0;
	end
    end
end