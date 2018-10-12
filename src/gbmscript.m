%function gbmscript(origRecon2)

%noiseCoeffs = 0:1;
%noiseCoeffs = [0];
%for i=1:length(noiseCoeffs)
%    [corefluxesmap,edgefluxesmap,coreexprmap,edgeexprmap,coretrueexprmap,edgetrueexprmap] = gbmfunc(noiseCoeffs(i),origRecon2);
%end

%function [corefluxesmap edgefluxesmap coreexprmap edgeexprmap coretrueexprmap edgetrueexprmap] = gbmfunc(noiseCoeff,origRecon2)

noiseCoeff = 0;
useGTEX = 1;
if 1
corefluxes = {};
edgefluxes = {};
coreids = {};
edgeids = {};
if useGTEX
outputDir1 = '/mnt/vdb/home/ubuntu2/MATLAB/FALCON/output/gbmscript/gtex/';
shortdata = textscan(fopen('/mnt/vdb/home/ubuntu2/MATLAB/FALCON/input/gtex/shortattributes.txt'),'%s%s%s%s%s%s','Delimiter','\t','HeaderLines',1);
phenodata = textscan(fopen('/mnt/vdb/home/ubuntu2/MATLAB/FALCON/input/gtex/GTEx_v7_Annotations_SubjectPhenotypesDS.txt'),'%s%s%s%s','Delimiter','\t','HeaderLines',1);
coreids = {};%'GTEX-12WSF','GTEX-13IVO','GTEX-13N1W','GTEX-13OVH','GTEX-13SLW','GTEX-148VJ','GTEX-14E7W','GTEX-14LLW','GTEX-14PK6','GTEX-15DCD'};
edgeids = {};%'GTEX-11EM3','GTEX-11NSD','GTEX-11P82','GTEX-11TT1','GTEX-12126','GTEX-12C56','GTEX-12WSE','GTEX-133LE','GTEX-13CZU','GTEX-13FTX'};
for i=1:length(shortdata{1})
    sampleshort = shortdata{1}{i};
    sampleshort = sampleshort(1:10);
    matchidx = find(strcmp(sampleshort,phenodata{1}));
    if ~isempty(matchidx) && strcmp(phenodata{3}{matchidx},'20-29') && strcmp(shortdata{3}{i},'Whole Blood') && strcmp(shortdata{5}{i},'TruSeq.v1') && strcmp(shortdata{6}{i},'RNASEQ') && isempty(regexp(shortdata{4}{i},'miR')) && (~isempty(regexp(shortdata{4}{i},'RNA isolation')) || ~isempty(regexp(shortdata{4}{i},'RNA Extraction')))
        coreids{end+1} = shortdata{1}{i};
    end
    matchidx = find(strcmp(sampleshort,phenodata{1}));
    if ~isempty(matchidx) && strcmp(phenodata{3}{matchidx},'70-79') && strcmp(shortdata{3}{i},'Whole Blood') && strcmp(shortdata{5}{i},'TruSeq.v1') && strcmp(shortdata{6}{i},'RNASEQ') && isempty(regexp(shortdata{4}{i},'miR')) && (~isempty(regexp(shortdata{4}{i},'RNA isolation')) || ~isempty(regexp(shortdata{4}{i},'RNA Extraction')))
        edgeids{end+1} = shortdata{1}{i};
    end
end
end
end

if 1
coreexpr = {};
edgeexpr = {};
runFlux = 1;
filenames = dir('/mnt/vdb/home/ubuntu2/MATLAB/FALCON/input/rnaseq1051');
if useGTEX
    filenames = dir('/mnt/vdb/home/ubuntu2/MATLAB/FALCON/input/gtex');
end
corefluxesmap = containers.Map;
edgefluxesmap = containers.Map;
coreexprmap = containers.Map;
edgeexprmap = containers.Map;
coretrueexprmap = containers.Map;
edgetrueexprmap = containers.Map;
for i=1:length(filenames)
    filename = filenames(i).name;
disp(i)
    if ~isempty(regexp(filename,'txt'))
        if useGTEX
	    matchesid = 0;
            matchescore = 0;
            for k=1:length(coreids)
		    if ~isempty(regexp(filename,coreids{k})) && length(keys(corefluxesmap))<10
	            matchesid = k;
                    matchescore = 1;
                end
	    end
	    for k=1:length(edgeids)
		if ~isempty(regexp(filename,edgeids{k})) && length(keys(edgefluxesmap))<10
	            matchesid = k;
                end
	    end
        end
        if useGTEX
	    if matchesid~=0
	        datafields = textscan(fopen(['/mnt/vdb/home/ubuntu2/MATLAB/FALCON/input/gtex/' filename]),'%s%f%f','HeaderLines',1);
            end
        else
            datafields = textscan(fopen(['/mnt/vdb/home/ubuntu2/MATLAB/FALCON/input/rnaseq1051/' filename]),'%s%f%f','HeaderLines',1);
        end
        if runFlux
	    if (useGTEX && matchesid~=0) || ~useGTEX
					          expData = datafields{2}/sum(datafields{2});
                trueExpData = expData;
                if noiseCoeff~=0
                    expData = expData.*abs(randn(length(expData),1)*noiseCoeff);
                    expData = expData/sum(expData)
                end
	        fluxes = runFluxMethod(expData,datafields{1},'test',constrainMediumExc(initializeRecon2(origRecon2)),'FALCON',datafields{3});
	        expr2 = obtainFALCONExp(expData,datafields{1},constrainMediumExc(initializeRecon2(origRecon2)),datafields{3});
	        trueexpr2 = obtainFALCONExp(trueExpData,datafields{1},constrainMediumExc(initializeRecon2(origRecon2)),datafields{3});
                if matchescore
                    corefluxesmap(coreids{matchesid}) = fluxes;
                    coreexprmap(coreids{matchesid}) = expr2;
                    coretrueexprmap(coreids{matchesid}) = trueexpr2;
                    outputDir2 = [outputDir1 'young/'];
                    corekey = coreids{matchesid};
	            writeData({fluxes,expr2,trueexpr2},[outputDir2 corekey '_noise' num2str(noiseCoeff)],'\t',{'flux','expr','trueexpr'});
                else
                    edgefluxesmap(edgeids{matchesid}) = fluxes;
                    edgeexprmap(edgeids{matchesid}) = expr2;
                    edgetrueexprmap(edgeids{matchesid}) = trueexpr2;
                    outputDir2 = [outputDir1 'old/'];
                    edgekey = edgeids{matchesid};
	            writeData({fluxes,expr2,trueexpr2},[outputDir2 corekey '_noise' num2str(noiseCoeff)],'\t',{'flux','expr','trueexpr'});
                end
		if length(keys(edgefluxesmap))>=10 && length(keys(corefluxesmap))>=10
		    nothing = 1;
		    %break
		end
            end
	else
	    expr = obtainFALCONExp(datafields{2},datafields{1},constrainMediumExc(initializeRecon2(origRecon2)),datafields{3});
            if ~isempty(regexp(filename,'core'))
                coreexpr{end+1} = expr;
            end
	    if ~isempty(regexp(filename,'edge'))
                edgeexpr{end+1} = expr;
            end
        end
    end
    if runFlux
        if ~useGTEX
	    if ~isempty(regexp(filename,'core'))
		coreids{end+1} = filename;
		corefluxes{end+1} = fluxes;
	    end
	    if ~isempty(regexp(filename,'edge'))
		edgeids{end+1} = filename;
		edgefluxes{end+1} = fluxes;
	    end
	end
    end
end
end

% exprpvalsarr = -ones(length(origRecon2.rxns),1);
% for i=1:1:length(origRecon2.rxns)
%     i
%     if all(~isnan(cellfun(@(x) x(i),edgeexpr))) && all(~isnan(cellfun(@(x) x(i),coreexpr)))
% 	pval = ranksum(cellfun(@(x) x(i),coreexpr),cellfun(@(x) x(i),edgeexpr));
% 	exprpvalsarr(i) = pval;
%     end
% end
