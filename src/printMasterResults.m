dataTable1 = {'xvals'};
dataTable2 = {'yvals'};
dataTable3 = {'xlabels'};

for i=1:length(reducedModelArr)
    dataTable1{end+1} = num2str(i);
    if isinf(ratioArr(i))
        dataTable2{end+1} = 'Inf';
    elseif isnan(ratioArr(i))
        dataTable2{end+1} = 'Not a Number';
    else
        dataTable2{end+1} = ratioArr(i);
    end
    if i==1
        dataTable3{end+1} = [num2str(i) '-' 'Base'];
    else
        dataTable3{end+1} = [num2str(i) '-' subsystemsAddedArr{i}];
    end
end

writeData({dataTable1,dataTable2,dataTable3},'/home/fs01/yw595/FALCONSeqGlucoseAll.txt','\t');