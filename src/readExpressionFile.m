function [ids data stds] = readExpressionFile(inputFile)

totalData = importdata(inputFile);
ids = totalData.data(1:end,1);
data = totalData.data(1:end,2);
stds = totalData.data(1:end,3);
ids = ids(~isnan(data));
data = data(~isnan(data));
stds = stds(~isnan(data));

end