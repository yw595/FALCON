if isunix
    FALCONDir = '/home/ubuntu/MATLAB/FALCON';
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
for i=1:length(origRecon2.genes)
    if ~isempty(regexp(origRecon2.genes{i},'\.'))
        %x = origRecon2.genes{i};
        origRecon2.genes{i} = strrep(origRecon2.genes{i},'.','');%x(1:regexp(x,'\.')-1);
    end
end
for i=1:length(origRecon2.grRules)
    if ~isempty(regexp(origRecon2.grRules{i},'\.'))
        %x = origRecon2.genes{i};
        origRecon2.grRules{i} = strrep(origRecon2.grRules{i},'.','');%x(1:regexp(x,'\.')-1);
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