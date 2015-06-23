load SigurdssonMouse.mat;
origSigurdssonMouse = mouse1415;

origLeeYeast = readCbModel('LeeYeast.xml');

load Recon2.v03.mat;
origRecon2 = modelRecon2beta121114_fixed_updatedID;
origRecon2.description = 'origRecon2';

origSigurdssonMouse.rxnNames = {};
for i=1:length(origSigurdssonMouse.rxns)
    matchHuman = strcmp( origRecon2.rxns,origSigurdssonMouse.rxns{i} );
    if any(matchHuman)
        origSigurdssonMouse.rxnNames{i} = origRecon2.rxnNames{matchHuman};
    else
        origSigurdssonMouse.rxnNames{i} = '';
    end
end
origSigurdssonMouse.rxnNames = origSigurdssonMouse.rxnNames';