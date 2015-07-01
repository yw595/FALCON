function statsArray = analyzeFlux(cellLineFile, cellLine, model)
% 
% INPUTS
%       cellLineFile - name of .mat file that contains predicted fluxes from a constraint-based method
%       cellLine - name of cell line in experimental CORE data, 
%       used its uptake/release profile calculate sensitivity wrt to
%       model - metabolic reconstruction, passed as argument to extractExcFlux
%
% OUTPUTS
%       statsArray - matrix of statistics measuring predicted versus experimental CORE flux,
%       each row corresponds to a different set of fluxes to measure,
%       each column to different statistic (see below for details)
%
% Author: Yiping Wang, 2015

[cellLinesArray, ~, coreTable, FVAVminArray, FVAVmaxArray] = readJainTable();

% get column in coreTable corresponding to cellLine, sort it by absolute flux size
coreTableCol = coreTable(:, strcmp(cellLinesArray, cellLine));
[sortedCoreTableCol, sortedCoreTableColIdxs] = sort(abs(coreTableCol), 1, 'descend');
sortedCoreTableCol = coreTableCol(sortedCoreTableColIdxs);

% based on cellLineFile, predicted flux vector has different names
% select right name, then extract predicted CORE exchange flux in v_Exc
% sort v_Exc, as well as FVAMin and Max for Supp Table 3
load(cellLineFile);
if regexp(cellLineFile,'fba')
    v_Exc = extractExcFlux(model, v_fba);
elseif regexp(cellLineFile,'gxfba')
    v_Exc = extractExcFlux(model, v_fba);
elseif regexp(cellLineFile,'eflux')
    v_Exc = extractExcFlux(model, v_fba);
else
    v_Exc = extractExcFlux(model, v_falcon);
end
sortedV_Exc = v_Exc(sortedCoreTableColIdxs);
sortedFVAVmaxArray = FVAVmaxArray(sortedCoreTableColIdxs);
sortedFVAVminArray = FVAVminArray(sortedCoreTableColIdxs);


%[columnVector(sortedCoreTableCol) columnVector(sortedCoreTableCol > 0) columnVector(sortedCoreTableCol < 0) columnVector(sortedFVAVmaxArray) columnVector(sortedFVAVminArray) columnVector(sortedV_Exc)]

% each row of statsArray represents a larger subset of CORE fluxes
% specifically, since the inner loop goes as j=1:i, the first row of statsArray represents
% stats for just the first, and largest, flux in sortedCoreTableCol
% the second row represents the two largest fluxes, and so on
% each column of statsArray contains one of 13 different statistics
% measuring divergence btw experimental and predicted CORE flux
statsArray = zeros( length(sortedCoreTableCol)+3,13 );
for i = 1:length(sortedCoreTableCol)
    uptakeTruePos=0;
    uptakeFalseNeg=0;
    releaseTruePos=0;
    releaseFalseNeg=0;
    includedIdxs=[];
    
    for j=1:i
        if sortedCoreTableCol(j) > 0
	    % check if model is capable of release for jth CORE flux
	    % through FVAVmax from Supp Table 3
	    % if yes, increment either releaseTruePos or FalseNeg
	    % similarly for uptake
            if sortedFVAVmaxArray(j) == 0
                statsArray(i,1) = statsArray(i,1)+1; %cannot match this CORE release based on FVA results
            else
                statsArray(i,2) = statsArray(i,2)+1; % can match this CORE release
                includedIdxs(end+1) = j;
                if(sortedV_Exc(j) > 0)
                    releaseTruePos = releaseTruePos + 1;
                else
                    releaseFalseNeg = releaseFalseNeg + 1;
                end
            end
        elseif sortedCoreTableCol(j) < 0
            if sortedFVAVminArray(j) == 0
                statsArray(i,3) = statsArray(i,3)+1; %cannot match this CORE uptake based on FVA results
            else
                statsArray(i,4) = statsArray(i,4)+1; %can match this CORE uptake
                includedIdxs(end+1) = j;
                if(sortedV_Exc(j) < 0)
                    uptakeTruePos = uptakeTruePos + 1;
                else
                    uptakeFalseNeg = uptakeFalseNeg + 1;
                end
            end
        else
            statsArray(i,5) = statsArray(i,5) + 1; % in case CORE flux is measured zero, therefore cannot assign true or false
            includedIdxs(end + 1) = j;
        end
    end

    % if uptake is happening when it shouldn't be (uptakeFalsePos),
    % then similarly release must not be happening when it should be, and similarly for other three
    uptakeFalsePos = releaseFalseNeg;
    uptakeTrueNeg = releaseTruePos;
    releaseFalsePos = uptakeFalseNeg;
    releaseTrueNeg = uptakeTruePos;
    
    if ~isempty(includedIdxs)
        Vest = columnVector(sortedV_Exc(includedIdxs));
        Vexp = columnVector(sortedCoreTableCol(includedIdxs));
        statsArray(i,6) = corr(Vest, Vexp, 'type', 'Pearson'); %Pearson correlation
        statsArray(i,7) = corr(Vest, Vexp, 'type', 'Spearman'); %Spearman correlation
        statsArray(i,8) = corr(Vest, Vexp, 'type', 'Kendall'); %Kendall correlation
        statsArray(i,9) = Vest' * Vexp / (norm(Vest) * norm(Vexp)); %cosine similarity
        statsArray(i,10) = norm(Vest - Vexp, 1); %L1 norm
        statsArray(i,11) = (uptakeTruePos + releaseTruePos) / ...
            (uptakeTruePos + uptakeFalseNeg + ...
            releaseTruePos + releaseFalseNeg); %overall sensitivity
        statsArray(i,12) = uptakeTruePos / ...
            (uptakeTruePos + uptakeFalseNeg); %uptake sensitivity
        statsArray(i,13) = releaseTruePos / ...
            (releaseTruePos + releaseFalseNeg); %release sensitivity
    end
end

summaryHandles = {@min, @max, @mean};
for i=1:size(statsArray,2)
    %filter for nans
    filteredStatsCol = statsArray( 1:length(sortedCoreTableCol),i );
    filteredStatsCol = filteredStatsCol(~isnan(filteredStatsCol));
    % append min, max and average rows for all stats
    % that is, for each statistic, output its min, max and mean value across all previous subsets,
    % like just the first flux, just first and second, and so forth
    for j=1:3
        if ~isempty(filteredStatsCol)
            statsArray( length(sortedCoreTableCol)+j,i ) = summaryHandles{j}(filteredStatsCol);
        else
            statsArray( length(sortedCoreTableCol)+j,i ) = nan;
        end
    end
end

end