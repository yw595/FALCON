if isunix
    FALCONDir = '/home/fs01/yw595/MATLAB/FALCON';
else
    FAlconDir = 'C:\Users\Yiping Wang\Documents\MATLAB\FALCON';
end
load([FALCONDir filesep 'input' filesep 'SigurdssonMouse.mat']);
origSigurdssonMouse = mouse1415;
origSigurdssonMouse.rxnNames = {};

origLeeYeast = readCbModel([FALCONDir filesep 'input' filesep 'LeeYeast.xml']);

load([FALCONDir filesep 'input' filesep 'Recon2.v03.mat']);
origRecon2 = modelRecon2beta121114_fixed_updatedID;
origRecon2.description = 'origRecon2';
fuckingStupid = 1;
for i=1:length(origRecon2.genes)
    if ~isempty(regexp(origRecon2.genes{i},'\.'))
        if fuckingStupid
            x = origRecon2.genes{i};
            origRecon2.genes{i} = x(1:regexp(x,'\.')-1);
        else
            origRecon2.genes{i} = strrep(origRecon2.genes{i},'.','');
        end
    end
end
for i=1:length(origRecon2.grRules)
    if ~isempty(regexp(origRecon2.grRules{i},'\.'))
        if fuckingStupid
            x = origRecon2.grRules{i};
            origRecon2.grRules{i} = x(1:regexp(x,'\.')-1);
        else
            origRecon2.grRules{i} = strrep(origRecon2.grRules{i},'.','');
        end
    end
end

for i=1:length(origSigurdssonMouse.rxns)
    matchHuman = strcmp( origRecon2.rxns,origSigurdssonMouse.rxns{i} );
    if any(matchHuman)
        origSigurdssonMouse.rxnNames{i} = origRecon2.rxnNames{matchHuman};
    else
        origSigurdssonMouse.rxnNames{i} = '';
    end
end
origSigurdssonMouse.rxnNames = origSigurdssonMouse.rxnNames;