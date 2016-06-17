function curatedRxnNames = makeCuratedRxnNames(model)

curatedGlycolysisRxnNames={'HEX1','PGI','PFK','FBA','TPI','GAPD','PGK','PGM','ENO','PYK','LDH_L'};
curatedPentosePhosphateRxnNames={'G6PDH2r','PGL','GND','RPI','RPE','TALA','TKT1','TKT2','PRPPS'};
curatedCitricAcidCycleRxnNames={'PDHm','CSm','ACONTm','ICDHxm','AKGDm','SUCOASm','SUCD1m','FUMm','MDHm'};
curatedOxidativePhosphorylationRxnNames={'NADH2_u10m','FADH2ETC','CYOR_u10m','CYOOm2','CYOOm3','ATPS4m','PPAm','PPA'};
curatedPurineSynthesisRxnNames={'GLUPRT','PRAGSr','GARFT','PRFGS','r0666','AIRCr','PRASCS','ADSL2','AICART','IMPC'};
curatedTransportMitoRxnNames={'PYRt2m','NADtm','FADtm','CO2tm','O2tm','H2Otm','PIt2m','ATPtm','O2Stm'};
curateddGTPSynthesisRxnNames={'IMPD','GMPS2','GK1','RNDR2','NDPK5'};
curateddATPSynthesisRxnNames={'ADSS','ADSL1','ADK1','RNDR1','NDPK8'};
curateddCTPSynthesisRxnNames={'CBPS','ASPCTr','DHORTS','DHORD9','ORPT','OMPDC','UMPK','NDPK2','CTPS2','NDPK3','RNDR3','NDPK7'};
curateddTTPSynthesisRxnNames={'CYTK2','DCMPDA','TMDS','DTMPK','NDPK4'};
curatedFolateMetabolismRxnNames={'DHFR','FTHFDH','MTHFC','MTHFD'};
curatedMiscRxnNames={'TRDR'};

curatedRxnNames=[curatedGlycolysisRxnNames';curatedPentosePhosphateRxnNames';curatedCitricAcidCycleRxnNames';curatedOxidativePhosphorylationRxnNames';curatedPurineSynthesisRxnNames';curatedTransportMitoRxnNames';curateddGTPSynthesisRxnNames';curateddATPSynthesisRxnNames';curateddCTPSynthesisRxnNames';curateddTTPSynthesisRxnNames';curatedFolateMetabolismRxnNames';curatedMiscRxnNames'];