function equation = printRxnEq(model,numberOrName)
    if isinteger(numberOrName)
        rxnIdx = numberOrName;
    else
        rxnIdx = find(strcmp(model.rxns,numberOrName));
    end
    
    inputMetIdxs = find(model.S(:,rxnIdx) < 0);
    outputMetIdxs = find(model.S(:,rxnIdx) > 0);
    equation = '';
    for i=1:length(inputMetIdxs)
        if i==1
	    equation = [equation '- '];
	end
        if i==length(inputMetIdxs)
	    equation = [equation num2str(abs(model.S(inputMetIdxs(i),rxnIdx))) ' ' model.mets{inputMetIdxs(i)} ' = '];
	else
	    equation = [equation num2str(abs(model.S(inputMetIdxs(i),rxnIdx))) ' ' model.mets{inputMetIdxs(i)} ' - '];
	end
    end

    for i=1:length(outputMetIdxs)
        if i==length(outputMetIdxs)
	    equation = [equation num2str(model.S(outputMetIdxs(i),rxnIdx)) ' ' model.mets{outputMetIdxs(i)}];
	else
	    equation = [equation num2str(model.S(outputMetIdxs(i),rxnIdx)) ' ' model.mets{outputMetIdxs(i)} ' + '];
	end
    end
end