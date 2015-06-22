inputDir = ['NCI60Sims' filesep 'nci60mRNA'];
inputFiles = dir(inputDir);

if ~exist([inputDir filesep 'outputReal'])
    system(['mkdir ' inputDir filesep 'outputReal']);
end

for i=1:length(inputFiles)
    %if ~strcmp(inputFiles(i).name, '.') && ~strcmp(inputFiles(i).name, '..') && ~strcmp(inputFiles(i).name, '.git') && ~strcmp(inputFiles(i).name, '.gitignore')
    if ~isempty(regexp(inputFiles(i).name, '.csv$'))
       runFALCONStripped(constrainMediumExc(initializeRecon2(origRecon2)), [inputDir filesep inputFiles(i).name], 1, [inputDir filesep 'output']);
    end
end

inputDir = ['NCI60Sims' filesep 'nci60prot'];
inputFiles = dir(inputDir);

if ~exist([inputDir filesep 'outputReal'])
    system(['mkdir ' inputDir filesep 'outputReal']);
end

for i=1:length(inputFiles)
    %if ~strcmp(inputFiles(i).name, '.') && ~strcmp(inputFiles(i).name, '..') && ~strcmp(inputFiles(i).name, '.git') && ~strcmp(inputFiles(i).name, '.gitignore')
    if ~isempty(regexp(inputFiles(i).name, '.csv$'))
       runFALCONStripped(constrainMediumExc(initializeRecon2(origRecon2)), [inputDir filesep inputFiles(i).name], 1, [inputDir filesep 'output']);
    end
end

inputDir = ['NCI60Sims' filesep 'nci60prot_mRNA'];
inputFiles = dir(inputDir);

if ~exist([inputDir filesep 'outputReal'])
    system(['mkdir ' inputDir filesep 'outputReal']);
end

for i=1:length(inputFiles)
    %if ~strcmp(inputFiles(i).name, '.') && ~strcmp(inputFiles(i).name, '..') && ~strcmp(inputFiles(i).name, '.git') && ~strcmp(inputFiles(i).name, '.gitignore')
    if ~isempty(regexp(inputFiles(i).name, '.csv$'))
       runFALCONStripped(constrainMediumExc(initializeRecon2(origRecon2)), [inputDir filesep inputFiles(i).name], 1, [inputDir filesep 'output']);
    end
end