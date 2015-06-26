function [cellLinesArray metsArray coreTable FVAVminArray FVAVmaxArray] = ...
    readJainTable(nomean)

warning('off', 'MATLAB:xlsread:ActiveX');
[excnumarray exctextarray raw] = xlsread( ... 
    ['input' filesep 'Supp Table 3 A community-driven global reconstruction ' ... 
     'of human metabolism 95.xls']);

cellLinesArray = exctextarray(9, 10:2:128);
metsArray = exctextarray(10:100, 2);
FVAVminArray = excnumarray(8:98, 1);
FVAVmaxArray = excnumarray(8:98, 3);

coreTableTemp = excnumarray(8:98, 8:127);

coreTable = [];
for i = 1:length(cellLinesArray)
    if exist('nomean','var') && nomean
        coreTable(:, (end+1):(end+2)) = coreTableTemp(:, (i*2-1):i*2);
    else
        coreTable(:, end + 1) = mean(coreTableTemp(:, (i*2-1):i*2), 2);
    end
end

end