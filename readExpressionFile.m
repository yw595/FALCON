function [ids data stds] = readExpressionFile(inputFile)

totalData = importdata(inputFile);
ids = totalData.data(2:end,1);
data = totalData.data(2:end,2);
stds = totalData.data(2:end,3);
ids = ids(~isnan(data));
data = data(~isnan(data));
stds = stds(~isnan(data));

end