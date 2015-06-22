function [cellLinesarray metsArray coreTable FVAVminArray FVAVmaxArray] = ...
    readJainTable(nomean)
%INPUT
% Assumes Jain's table is in the running directory
%
%OUTPUT
%
%
warning('off', 'MATLAB:xlsread:ActiveX');
[excnumarray exctextarray raw] = xlsread( ... 
    ['Supp Table 3 A community-driven global reconstruction ' ... 
     'of human metabolism 95.xls']);
[height width] = size(excnumarray);
coreTable1 = excnumarray(8:98, 8:width);
cellLinesArray1 = exctextarray(9, 10:2:128);
coreTable = [];
cellLinesArray = {};

for i = 1:length(cellLinesArray1)
    cellLinesArray{end + 1} = cellLinesArray1{i};
    if exist('nomean','var')
        if nomean
            % In case we want to look at individual replicates
            coreTable(:, (end+1):(end+2)) = coreTable1(:, (i*2-1):i*2);
        else
            coreTable(:, end + 1) = mean(coreTable1(:, (i*2-1):i*2), 2);
        end
    else
      coreTable(:, end + 1) = mean(coreTable1(:, (i*2-1):i*2), 2);
    end
end
metsArray = exctextarray(10:100, 2);
FVAVminArray = excnumarray(8:98, 1);
FVAVmaxArray = excnumarray(8:98, 3);
end

