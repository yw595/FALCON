function modelToRun = analyzeFluxScriptHelper(cellLine, inputPrefix)
    load(['NCI60Sims' filesep 'nci60prot' filesep 'specificModelsPar' ...
                    filesep 'specificModel' cellLine inputPrefix '.mat']);
                    eval(['modelToRun = specificModel' inputPrefix ';']);
end