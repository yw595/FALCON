nanArr = []; run = [];
for z=1:length(ratioArr2Arr)
    ratioArr2 = ratioArr2Arr{z};
    nanArr = [nanArr isnan(ratioArr2) | isinf(ratioArr2)];
    run = [run z*ones(1,length(ratioArr2))];
end

writeData({nanArr,run},'/home/fs01/yw595/nanInfCDFs.txt','\t',{'nan','run'});

single = 0;
if single
    ratioArr2Arr = {ratioArr2};
end
nanInfPercentsArr = [];

for z=1:length(ratioArr2Arr)
ratioArr2 = ratioArr2Arr{z};
nanInfBins = histc(find(isnan(ratioArr2) | isinf(ratioArr2)),[0,100,200,300,400,500,600,700,Inf]);
nanInfLabels = {}; nanInfPercents = [];
for i=1:length(nanInfBins)-2
    nanInfLabels{i} = [num2str(100*(i-1)+1) '-' num2str(100*i)];
    nanInfPercents(i) = nanInfBins(i)/100;
end
nanInfLabels{end+1} = [num2str(100*(length(nanInfLabels))+1) '-' num2str(length(ratioArr2))];
nanInfPercents(end+1) = nanInfBins(end-1)/mod(length(ratioArr2),100);
nanInfPercentsArr(:,z) = nanInfPercents;
end

if single
nanInfPercents = nanInfPercentsArr(:,1);
writeData({nanInfLabels,nanInfPercents},[transferDir '/nanInfBins.txt'],'\t',{'reactionNumber','percentage'});
else
    nanInfPercents = mean(nanInfPercentsArr,2);
    nanInfPercentsLower = min(nanInfPercentsArr,[],2);
    nanInfPercentsUpper = max(nanInfPercentsArr,[],2);
    writeData({nanInfLabels,nanInfPercents,nanInfPercentsLower,nanInfPercentsUpper},['/home/fs01/yw595/nanInfBinsErrBar.txt'],'\t',{'reactionNumber','percentage','lower','upper'});
end



