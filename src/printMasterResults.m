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

dataTable1 = {'xvals'};
dataTable2 = {'yvals'};
dataTable3 = {'xlabels'};

for i=1:length(reducedModelArrSkip)
    dataTable1{end+1} = num2str(i);
    if isinf(ratioArrSkip(i))
        dataTable2{end+1} = 'Inf';
    elseif isnan(ratioArrSkip(i))
        dataTable2{end+1} = 'Not a Number';
    else
        dataTable2{end+1} = ratioArrSkip(i);
    end
    if i==1
        dataTable3{end+1} = [num2str(i) '-' 'Base'];
    else
        dataTable3{end+1} = [num2str(i) '-' subsystemsAddedArrSkip{i}];
    end
end

writeData({dataTable1,dataTable2,dataTable3},'/home/fs01/yw595/FALCONSeqGlucoseAllSkip.txt','\t');

for k=1:ceil(length(ratioArr2)/100)
    dataTable1 = {'xvals'};
    dataTable2 = {'yvals'};
    dataTable3 = {'xlabels'};

    for i=(k-1)*100+1:min(length(ratioArr2),k*100)
        dataTable1{end+1} = num2str(i);
        if isinf(ratioArr2(i))
            dataTable2{end+1} = 'Inf';
        elseif isnan(ratioArr2(i))
            dataTable2{end+1} = 'Not a Number';
        else
            dataTable2{end+1} = ratioArr2(i);
        end
        if i==1
            dataTable3{end+1} = [num2str(i) '-' 'Base'];
        else
            dataTable3{end+1} = [num2str(i) '-' reactionsToAdd{i}];
        end
    end

    writeData({dataTable1,dataTable2,dataTable3},['/home/fs01/yw595/FALCONSeqGlucoseAllExchange' num2str(k) '.txt'],'\t');

end
