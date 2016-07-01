inputFiles = dir('HDSims/HDFalconData');

for i=1:length(inputFiles)
    if ~strcmp(inputFiles(i).name, '.') && ~strcmp(inputFiles(i).name, '..')
       runFALCONStripped(newSigurdssonMouse,['HDSims/HDFalconData/' inputFiles(i).name],1,'HDSims');
    end
end