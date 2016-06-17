function [transportMets demandMets biomassMets] = makeToAddMets()

% add NOTE: FUNCTIONAL TRANSPORt & EXCHANGE RXNS, EVEN IF TF_ NAME 
% Rxns for transportMets, NOTE: DOES INCLUDE GLC_D[C] ALSO,
% IN EARLY CASES WHERE NO OTHER GLUCOSE YET ADDED, DOESN'T SHOW
% UP IN SIMPLEREDUCEDMODELSCRIPT FIGURE
% add biomass drain for demandMets, drain rxns for demandMets2
transportMets={'o2[c]','h2o[c]','h[c]','pi[c]','glc_D[c]','lac_L[c]','co2[c]','gln_L[c]','glu_L[c]','gly[c]','asp_L[c]','arg_L[c]','asn_L[c]','cys_L[c]','his_L[c]','4hpro_LT[c]','ile_L[c]','leu_L[c]','lys_L[c]','met_L[c]','phe_L[c]','pro_L[c]','ser_L[c]','thr_L[c]','trp_L[c]','tyr_L[c]','val_L[c]','btn[c]','chol[c]','pnto_R[c]','fol[c]','ncam[c]','bz[c]','pydxn[c]','ribflv[c]','thm[c]','adpcbl[c]','inost[c]','ca2[c]','so4[c]','k[e]','cl[c]','na1[c]','gthrd[c]','nad[c]','nadh[c]','fum[c]','o2s[c]','hco3[c]','q10[m]','q10h2[m]','nh4[c]','gtp[c]','gdp[c]'};
demandMets={'gln_L[c]','glu_L[c]','gly[c]','asp_L[c]','arg_L[c]', 'asn_L[c]','cys_L[c]','his_L[c]','4hpro_LT[c]','ile_L[c]','leu_L[c]','lys_L[c]','met_L[c]','phe_L[c]','pro_L[c]','ser_L[c]','thr_L[c]','trp_L[c]','tyr_L[c]','val_L[c]'};
biomassMets={'dctp[c]','dttp[c]','dgtp[c]','datp[c]'};