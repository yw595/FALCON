function appendArr = appendNumFunc(arr)

maxLen = length(num2str(length(arr)));
appendArr = arr;
for i=1:length(appendArr)
    appendNum = num2str(i);
    while length(appendNum)<maxLen
        appendNum = ['0' appendNum];
    end
    appendArr{i} = [appendNum '_' appendArr{i}];
end

end

