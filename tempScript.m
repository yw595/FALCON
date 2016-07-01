%makeTissueSpecificModelsNCI;
%nci60Script;
matlabpool close;
matlabpool 8;
nci60Script;
analyzeFluxScript;
averageScript;