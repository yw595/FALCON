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
tissuestomeasure = unique(values(idstotissue));

for z1=7:length(tissuestomeasure)
    coreidszth = {};
    coreorigexprzth = {};
    for i=1:length(coreids)
	if strcmp(idstotissue(coreids{i}),tissuestomeasure{z1})
	    coreidszth{end+1} = coreids{i};
            coreorigexprzth{end+1} = coreorigexpr{i};
        end
    end
    edgeidszth = {};
    edgeorigexprzth = {};
    for i=1:length(edgeids)
	if strcmp(idstotissue(edgeids{i}),tissuestomeasure{z1})
	    edgeidszth{end+1} = edgeids{i};
            edgeorigexprzth{end+1} = edgeorigexpr{i};
        end
    end

    if length(edgeorigexprzth)==0 || length(coreorigexprzth)==0
	  continue;
    end
	
    coremean = [];
    for i=1:length(coreorigexprzth{1})
	ithExpr = [];
	for j=1:length(coreorigexprzth)
	    ithExpr(j) = coreorigexprzth{j}(i);
	end
	if ~any(isnan(ithExpr))
	    coremean(i) = mean(ithExpr);
	else
	    coremean(i) = NaN;
	end
    end
    edgemean = [];
    for i=1:length(edgeorigexprzth{1})
	ithExpr = [];
	for j=1:length(edgeorigexprzth)
	    ithExpr(j) = edgeorigexprzth{j}(i);
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
	tissuestomeasure{z1}
	z
	randvec = [];
	for i=1:length(edgeorigexprzth{1})
	    ithExpr = [];
	    if useedge
		for j=1:length(edgeorigexprzth)
		    ithExpr(end+1) = edgeorigexprzth{j}(i);
		end
	    end
	    if usecore
		for j=1:length(coreorigexprzth)
		    ithExpr(end+1) = coreorigexprzth{j}(i);
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
            writeData({randvec},[outputDir1 '/edgerandexpr' tissuestomeasure{z1} num2str(z)],'\t');
	end
	if useedge==0 && usecore==1
	    randvec = randvec-coremean;
            writeData({randvec},[outputDir1 '/corerandexpr' tissuestomeasure{z1} num2str(z)],'\t');
	end
	if useedge==1 && usecore==1
	    allmean = (edgemean*length(edgeorigexprzth)+coremean*length(coreorigexprzth))/(length(edgeorigexprzth)+length(coreorigexprzth));
	    randvec = randvec-allmean;
            writeData({randvec},[outputDir1 '/allrandexpr' tissuestomeasure{z1} num2str(z)],'\t');
	end
    end
end
