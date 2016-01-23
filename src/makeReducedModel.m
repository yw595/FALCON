function reducedModel=makeReducedModel(model,toAddRxnNames,description)
if ~exist('toAddRxnNames','var')
glycolysisRxnNames=model.rxns(strcmp(model.subSystems,'Glycolysis/gluconeogenesis'));
pentosePhosphateRxnNames=model.rxns(strcmp(model.subSystems,'Pentose phosphate pathway'));
purineSynthesisRxnNames=model.rxns(strcmp(model.subSystems,'Purine synthesis'));
pyrimidineSynthesisRxnNames=model.rxns(strcmp(model.subSystems,'Pyrimidine synthesis'));
citricAcidCycleRxnNames=model.rxns(strcmp(model.subSystems,'Citric acid cycle'));
oxidativePhosphorylationRxnNames=model.rxns(strcmp(model.subSystems,'Oxidative phosphorylation'));
exchangeDemandRxnNames=model.rxns(strcmp(model.subSystems,'Exchange/demand reaction'));
transportMitoRxnNames=model.rxns(strcmp(model.subSystems,'Transport, mitochondrial'));
transportExRxnNames=model.rxns(strcmp(model.subSystems,'Transport, extracellular'));
transportNuclearRxnNames=model.rxns(strcmp(model.subSystems,'Transport, nuclear'));
nucleotideInterconversionRxnNames=model.rxns(strcmp(model.subSystems,'Nucleotide interconversion'));
fattyAcidSynthesisRxnNames=model.rxns(strcmp(model.subSystems,'Fatty acid synthesis'));
purineCatabolismRxnNames=model.rxns(strcmp(model.subSystems,'Purine catabolism'));
pyrimidineCatabolismRxnNames=model.rxns(strcmp(model.subSystems,'Pyrimidine catabolism'));
folateMetabolismRxnNames=model.rxns(strcmp(model.subSystems,'Folate metabolism'));
excludeRxnNames={exchangeDemandRxnNames; transportExRxnNames; purineCatabolismRxnNames; pyrimidineCatabolismRxnNames};
toAddRxnNames={glycolysisRxnNames; pentosePhosphateRxnNames; purineSynthesisRxnNames;};% pyrimidineSynthesisRxnNames; oxidativePhosphorylationRxnNames; citricAcidCycleRxnNames; nucleotideInterconversionRxnNames; folateMetabolismRxnNames};
end

if (~exist('description','var'))
  description='';
end

% at bottom of block, curatedRxnNames adds together all
simpleGlycolysisRxnNames={'HEX1','PGI','PFK','FBA','TPI','GAPD','PGK','PGM','ENO','PYK','LDH_L'};
simplePentosePhosphateRxnNames={'G6PDH2r','PGL','GND','RPI','RPE','TALA','TKT1','TKT2','PRPPS'};
simpleCitricAcidCycleRxnNames={'PDHm','CSm','ACONTm','ICDHxm','AKGDm','SUCOASm','SUCD1m','FUMm','MDHm'};
simpleOxidativePhosphorylationRxnNames={'NADH2_u10m','FADH2ETC','CYOR_u10m','CYOOm2','CYOOm3','ATPS4m','PPAm','PPA'};
simplePurineSynthesisRxnNames={'GLUPRT','PRAGSr','GARFT','PRFGS','r0666','AIRCr','PRASCS','ADSL2','AICART','IMPC'};
simpleTransportMitoRxnNames={'PYRt2m','NADtm','FADtm','CO2tm','O2tm','H2Otm','PIt2m','ATPtm','O2Stm'};
simpledGTPSynthesisRxnNames={'IMPD','GMPS2','GK1','RNDR2','NDPK5'};
simpledATPSynthesisRxnNames={'ADSS','ADSL1','ADK1','RNDR1','NDPK8'};
simpledCTPSynthesisRxnNames={'CBPS','ASPCTr','DHORTS','DHORD9','ORPT','OMPDC','UMPK','NDPK2','CTPS2','NDPK3','RNDR3','NDPK7'};
simpledTTPSynthesisRxnNames={'CYTK2','DCMPDA','TMDS','DTMPK','NDPK4'};
simpleFolateMetabolismRxnNames={'DHFR','FTHFDH','MTHFC','MTHFD'};
simpleMiscRxnNames={'TRDR'};
curatedRxnNames=[simpleGlycolysisRxnNames';simplePentosePhosphateRxnNames';simpleCitricAcidCycleRxnNames';simpleOxidativePhosphorylationRxnNames';simplePurineSynthesisRxnNames';simpleTransportMitoRxnNames';simpledGTPSynthesisRxnNames';simpledATPSynthesisRxnNames';simpledCTPSynthesisRxnNames';simpledTTPSynthesisRxnNames';simpleFolateMetabolismRxnNames';simpleMiscRxnNames'];

% assume toAddRxnNames is cell array of cell arrays, add rxnName if not in curatedRxnNames
for i=1:length(toAddRxnNames)
  toAddRxnNamesSet=toAddRxnNames{i};
  for j=1:length(toAddRxnNamesSet)
    if(sum(strcmp(toAddRxnNamesSet{j},curatedRxnNames))==0)
      curatedRxnNames=[curatedRxnNames; {toAddRxnNamesSet{j}}];
    end
  end
end

% subselect rxns, mets and S
selRxns = ismember(model.rxns,curatedRxnNames);
subS = model.S(:,selRxns);
if (nargin < 8)
    selMets = ~all(subS == 0,2);
else
    selMets = ismember(model.mets,metNames);
end

subS = subS(selMets,:);

subModel.S = subS;
subModel.rxns = model.rxns(selRxns);
subModel.mets = model.mets(selMets);

% subselect b, metNames, Formulas, rev, lb, ub, c, genes (selGenes are any associated with selRxns),
% geneNames, rxnNames, subSystems
if (isfield(model,'b'))
    subModel.b = model.b(selMets);
end
if (isfield(model,'metNames'))
    subModel.metNames = model.metNames(selMets);
end
if (isfield(model,'metFormulas'))
    subModel.metFormulas = model.metFormulas(selMets);
end
if (isfield(model,'description'))
    subModel.description = description;
end
if (isfield(model,'rev'))
    subModel.rev = model.rev(selRxns);
end
if (isfield(model,'lb'))
    subModel.lb = model.lb(selRxns);
end 
if (isfield(model,'ub'))
    subModel.ub = model.ub(selRxns);
end
if (isfield(model,'c'))
    subModel.c = model.c(selRxns);
end
if (isfield(model,'genes'))
   newRxnGeneMat = model.rxnGeneMat(selRxns,:); 
   selGenes = sum(newRxnGeneMat)' > 0;
   subModel.rxnGeneMat = newRxnGeneMat(:,selGenes);
   subModel.genes = model.genes(selGenes);
   subModel.grRules = model.grRules(selRxns);
   subModel.rules = model.rules(selRxns);
end
if (isfield(model,'geneNames'))
    subModel.geneNameRules = model.geneNameRules(selRxns);
    subModel.geneNames = model.geneNames(selGenes);
end
if (isfield(model,'subSystems'))
    subModel.subSystems = model.subSystems(selRxns);
end
if(isfield(model,'rxnNames'))
    subModel.rxnNames=model.rxnNames(selRxns);
end
if isfield(model,'reversibleModel')
    subModel.reversibleModel=model.reversibleModel;
end

% add NOTE: FUNCTIONAL TRANSPORt & EXCHANGE RXNS, EVEN IF TF_ NAME 
% Rxns for transportMets, NOTE: DOES INCLUDE GLC_D[C] ALSO,
% IN EARLY CASES WHERE NO OTHER GLUCOSE YET ADDED, DOESN'T SHOW
% UP IN SIMPLEREDUCEDMODELSCRIPT FIGURE
% add biomass drain for demandMets, drain rxns for demandMets2
reducedModel=subModel;
transportMets={'o2[c]','h2o[c]','h[c]','pi[c]','glc_D[c]','lac_L[c]','co2[c]','gln_L[c]','glu_L[c]','gly[c]','asp_L[c]','arg_L[c]','asn_L[c]','cys_L[c]','his_L[c]','4hpro_LT[c]','ile_L[c]','leu_L[c]','lys_L[c]','met_L[c]','phe_L[c]','pro_L[c]','ser_L[c]','thr_L[c]','trp_L[c]','tyr_L[c]','val_L[c]','btn[c]','chol[c]','pnto_R[c]','fol[c]','ncam[c]','bz[c]','pydxn[c]','ribflv[c]','thm[c]','adpcbl[c]','inost[c]','ca2[c]','so4[c]','k[e]','cl[c]','na1[c]','gthrd[c]','nad[c]','nadh[c]','fum[c]','o2s[c]','hco3[c]','q10[m]','q10h2[m]','nh4[c]','gtp[c]','gdp[c]'};
demandMets={'dctp[c]','dttp[c]','dgtp[c]','datp[c]'};
demandMets2={'gln_L[c]','glu_L[c]','gly[c]','asp_L[c]','arg_L[c]','asn_L[c]','cys_L[c]','his_L[c]','4hpro_LT[c]','ile_L[c]','leu_L[c]','lys_L[c]','met_L[c]','phe_L[c]','pro_L[c]','ser_L[c]','thr_L[c]','trp_L[c]','tyr_L[c]','val_L[c]'};
transportRxns={}; 

% fix four reactions that have weird bounds in origRecon2, NOTE: ALSO PYK LOWER LIM SET TO 1
for i=1:length(reducedModel.rxns)
    if(strcmp(reducedModel.rxns{i},'PYK'))
        reducedModel.lb(i)=1;
	reducedModel.ub(i)=1000;
	reducedModel.rev(i)=0;
    end
    if(strcmp(reducedModel.rxns{i},'PIt2m'))
        reducedModel.lb(i)=-1000;
	reducedModel.ub(i)=1000;
	reducedModel.rev(i)=1;
    end
    if(strcmp(reducedModel.rxns{i},'CYOOm2'))
        reducedModel.lb(i)=-1000;
	reducedModel.ub(i)=1000;
	reducedModel.rev(i)=1;
    end
    if(strcmp(reducedModel.rxns{i},'FTHFDH'))
        reducedModel.lb(i)=-1000;
	reducedModel.ub(i)=1000;
	reducedModel.rev(i)=1;
    end
end

for i=1:length(reducedModel.mets)
    internalIdx=i;
    if(sum(strcmp([reducedModel.mets{i}],transportMets))~=0)
        reducedModel.rxns{end+1}=['TR_' reducedModel.mets{i}];
	reducedModel.S(internalIdx,end+1)=1;
	reducedModel.rev(end+1)=1;
	reducedModel.lb(end+1)=-1000;
	reducedModel.ub(end+1)=1000;
	reducedModel.c(end+1)=0;
	reducedModel.rxnGeneMat(end+1,:)=zeros(1,size(reducedModel.rxnGeneMat,2));
	reducedModel.grRules{end+1}='';
	reducedModel.subSystems{end+1}='';
	reducedModel.rxnNames{end+1}='';
	reducedModel.rules{end+1}='';
    end
    if(sum(strcmp([reducedModel.mets{i}],demandMets2))~=0)
        reducedModel.rxns{end+1}=['DEMAND_' reducedModel.mets{i}];
	reducedModel.S(internalIdx,end+1)=-1;
	reducedModel.rev(end+1)=0;
	reducedModel.lb(end+1)=0;
	reducedModel.ub(end+1)=1000;
	reducedModel.c(end+1)=0;
	reducedModel.rxnGeneMat(end+1,:)=zeros(1,size(reducedModel.rxnGeneMat,2));
	reducedModel.grRules{end+1}='';
	reducedModel.subSystems{end+1}='';
	reducedModel.rxnNames{end+1}='';
	reducedModel.rules{end+1}='';
    end
end

reducedModel.rxns{end+1}=['BIOMASS'];
reducedModel.rev(end+1)=0;
reducedModel.lb(end+1)=0;
reducedModel.ub(end+1)=1000;
reducedModel.c(end+1)=0;
reducedModel.rxnGeneMat(end+1,:)=zeros(1,size(reducedModel.rxnGeneMat,2));
reducedModel.grRules{end+1}='';
reducedModel.subSystems{end+1}='';
reducedModel.rxnNames{end+1}='';
reducedModel.rules{end+1}='';
added=0;
for i=1:length(reducedModel.mets)
if(sum(strcmp([reducedModel.mets{i}],demandMets))~=0)
    internalIdx=i;
    if(added==0)
        reducedModel.S(internalIdx,end+1)=-1;
	added=1;
    else
        reducedModel.S(internalIdx,end)=-1;
    end
end
end