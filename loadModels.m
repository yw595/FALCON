load(['input' filesep 'SigurdssonMouse.mat']);
origSigurdssonMouse = mouse1415;
origSigurdssonMouse.rxnNames = {};

origLeeYeast = readCbModel(['input' filesep 'LeeYeast.xml']);

load(['input' filesep 'Recon2.v03.mat']);
origRecon2 = modelRecon2beta121114_fixed_updatedID;
origRecon2.description = 'origRecon2';

for i=1:length(origSigurdssonMouse.rxns)
    matchHuman = strcmp( origRecon2.rxns,origSigurdssonMouse.rxns{i} );
    if any(matchHuman)
        origSigurdssonMouse.rxnNames{i} = origRecon2.rxnNames{matchHuman};
    else
        origSigurdssonMouse.rxnNames{i} = '';
    end
end
origSigurdssonMouse.rxnNames = origSigurdssonMouse.rxnNames;