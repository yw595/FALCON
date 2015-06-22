inputDir = ['NCI60Sims' filesep 'nci60mRNA'];
inputFiles = dir(inputDir);

for i=1:10
    if ~strcmp(inputFiles(i).name, '.') && ~strcmp(inputFiles(i).name, '..')
       inputFI = fopen([inputDir filesep inputFiles(i).name]);
       outputFI = fopen([inputDir filesep inputFiles(i).name '.temp'],'w');
       line = fgetl(inputFI);
       while line~=-1
       	   if isempty(regexp(line,'nan'))
       	       fprintf(outputFI,[line '\n']);
	   end
	   line = fgetl(inputFI);
       end
       fclose(outputFI);
       fclose(inputFI);
       system(['rm ' inputDir filesep inputFiles(i).name]);
       system(['mv ' inputDir filesep inputFiles(i).name '.temp' ' ' inputDir filesep inputFiles(i).name])
    end
end