function statsArray = analyzeFlux(cellLineFile, cellLine, model)

[cellLinesArray, ~, coreTable, FVAVminArray, FVAVmaxArray] = readJainTable();

coreTableCol = coreTable(:, strcmp(cellLinesArray, cellLine));
[sortedCoreTableCol, sortedCoreTableColIdxs] = sort(abs(coreTableCol), 1, 'descend');
load(cellLineFile);
sortedCoreTableCol = coreTableCol(sortedCoreTableColIdxs);
v_Exc = extractExcFlux(model, v_falcon);
sortedV_Exc = v_Exc(sortedCoreTableColIdxs);
sortedFVAVmaxArray = FVAVmaxArray(sortedCoreTableColIdxs);
sortedFVAVminArray = FVAVminArray(sortedCoreTableColIdxs);


[columnVector(sortedCoreTableCol) columnVector(sortedCoreTableCol > 0) columnVector(sortedCoreTableCol < 0) columnVector(sortedFVAVmaxArray) columnVector(sortedFVAVminArray) columnVector(sortedV_Exc)]

statsArray = zeros( length(sortedCoreTableCol)+3,13 );
for i = 1:length(sortedCoreTableCol)
    uptakeTruePos=0;
    uptakeFalseNeg=0;
    releaseTruePos=0;
    releaseFalseNeg=0;
    includedIdxs=[];
    
    for j=1:i
        if sortedCoreTableCol(j) > 0
            if sortedFVAVmaxArray(j) == 0
                statsArray(i,1) = statsArray(i,1)+1;
            else
                statsArray(i,2) = statsArray(i,2)+1;
                includedIdxs(end+1) = j;
                if(sortedV_Exc(j) > 0)
                    releaseTruePos = releaseTruePos + 1;
                else
                    releaseFalseNeg = releaseFalseNeg + 1;
                end
            end
        elseif sortedCoreTableCol(j) < 0
            if sortedFVAVminArray(j) == 0
                statsArray(i,3) = statsArray(i,3)+1;
            else
                statsArray(i,4) = statsArray(i,4)+1;
                includedIdxs(end+1) = j;
                if(sortedV_Exc(j) < 0)
                    uptakeTruePos = uptakeTruePos + 1;
                else
                    uptakeFalseNeg = uptakeFalseNeg + 1;
                end
            end
        else
            statsArray(i,5) = statsArray(i,5) + 1;
            includedIdxs(end + 1) = j;
        end
    end
    %releaseFalseNeg
    %releaseTruePos
    %uptakeFalseNeg
    %uptakeTruePos

    uptakeFalsePos = releaseFalseNeg;
    uptakeTrueNeg = releaseTruePos;
    releaseFalsePos = uptakeFalseNeg;
    releaseTrueNeg = uptakeTruePos;
    
    if ~isempty(includedIdxs)
        Vest = columnVector(sortedV_Exc(includedIdxs));
        Vexp = columnVector(sortedCoreTableCol(includedIdxs));
        statsArray(i,6) = corr(Vest, Vexp, 'type', 'Pearson');
        statsArray(i,7) = corr(Vest, Vexp, 'type', 'Spearman');
        statsArray(i,8) = corr(Vest, Vexp, 'type', 'Kendall');
        statsArray(i,9) = Vest' * Vexp / (norm(Vest) * norm(Vexp));
        statsArray(i,10) = norm(Vest - Vexp, 1);
        statsArray(i,11) = (uptakeTruePos + releaseTruePos) / ...
            (uptakeTruePos + uptakeFalseNeg + ...
            releaseTruePos + releaseFalseNeg);
        statsArray(i,12) = uptakeTruePos / ...
            (uptakeTruePos + uptakeFalseNeg);
        statsArray(i,13) = releaseTruePos / ...
            (releaseTruePos + releaseFalseNeg);
    end
end

summaryHandles = {@min, @max, @mean};
for i=1:size(statsArray,2)
    filteredStatsCol = statsArray( 1:length(sortedCoreTableCol),i );
    filteredStatsCol = filteredStatsCol(~isnan(filteredStatsCol));
    for j=1:3
        if ~isempty(filteredStatsCol)
            statsArray( length(sortedCoreTableCol)+j,i ) = summaryHandles{j}(filteredStatsCol);
        else
            statsArray( length(sortedCoreTableCol)+j,i ) = nan;
        end
    end
end

end