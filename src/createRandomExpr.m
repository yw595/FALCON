if 0
corrmat = -ones(7440,7440);
for i=1:7440
    i
    ithExpr = [];
    for k=1:length(coreexpr)
	ithExpr(k) = coreexpr{k}(i);
    end
    for j=1:7440
	jthExpr = [];
        for k=1:length(coreexpr)
	    jthExpr(k) = coreexpr{k}(j);
        end
        if ~any(isnan(ithExpr)) && ~any(isnan(jthExpr)) && min(ithExpr)~=max(ithExpr) && min(jthExpr)~=max(jthExpr)
            [rho, pval] = corr(ithExpr',jthExpr');
            rho
            corrmat(i,j) = rho*rho;
        end
    end
end
end

outputDir1 = '/mnt/vdb/home/ubuntu2/MATLAB/FALCON/output/createRandomExpr';
if ~exist(outputDir1,'dir')
    mkdir(outputDir1);
end
coremean = [];
for i=1:length(coreorigexpr{1})
    ithExpr = [];
    for j=1:length(coreorigexpr)
	ithExpr(j) = coreorigexpr{j}(i);
    end
    if ~any(isnan(ithExpr))
        coremean(i) = mean(ithExpr);
    else
        coremean(i) = NaN;
    end
end
edgemean = [];
for i=1:length(edgeorigexpr{1})
    ithExpr = [];
    for j=1:length(edgeorigexpr)
	ithExpr(j) = edgeorigexpr{j}(i);
    end
    if ~any(isnan(ithExpr))
        edgemean(i) = mean(ithExpr);
    else
        edgemean(i) = NaN;
    end
end
    usecore = 1;
useedge = 1;
for z=1:100
    z
    randvec = [];
    for i=1:length(edgeorigexpr{1})
	ithExpr = [];
        if useedge
	    for j=1:length(edgeorigexpr)
		ithExpr(end+1) = edgeorigexpr{j}(i);
	    end
	end
	if usecore
	    for j=1:length(coreorigexpr)
		ithExpr(end+1) = coreorigexpr{j}(i);
	    end
	end
	
	if ~any(isnan(ithExpr))
	    randvec(i) = quantile(ithExpr,rand(1));
	else
	    randvec(i) = NaN;
	end
    end
    if useedge==1 && usecore==0
        randvec = randvec-edgemean;
        writeData({randvec},[outputDir1 '/edgerandexpr' num2str(z)],'\t');
    end
    if useedge==0 && usecore==1
        randvec = randvec-coremean;
        writeData({randvec},[outputDir1 '/corerandexpr' num2str(z)],'\t');
    end
    if useedge==1 && usecore==1
      allmean = (edgemean*length(edgeorigexpr)+coremean*length(coreorigexpr))/(length(edgeorigexpr)+length(coreorigexpr));
        randvec = randvec-allmean;
        writeData({randvec},[outputDir1 '/allrandexpr' num2str(z)],'\t');
    end
end
