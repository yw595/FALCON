function [originTissuesArray INITFilesArray mCADREFilesArray ] = makeOriginTissuesArray(cellLinesArray)
    originTissuesArray = {}; INITFilesArray = {}; mCADREFilesArray = {};
    for i=1:length(cellLinesArray);
        if any(strcmp(cellLinesArray{i},{'NCI_H23','NCI_H522','A549_ATCC','EKVX','NCI_H226','NCI_H322M','NCI_H460','HOP_62','HOP_92'}))
	    originTissuesArray{end+1} = 'lung';
	    INITFilesArray{end+1} = 'iLungCancer1490';
	    mCADREFilesArray{end+1} = 'lung_tumor';
	elseif any(strcmp(cellLinesArray{i},{'HT29','HCC_2998','HCT_116','SW620','COLO_205','HCT_15','KM12'}))
	    originTissuesArray{end+1} = 'colon';
	    INITFilesArray{end+1} = 'iColorectalCancer1750';
	    mCADREFilesArray{end+1} = 'colon_tumor';
	elseif any(strcmp(cellLinesArray{i},{'MCF7','NCI_ADR_RES','MDA_MB_231_ATCC','HS_578T','MDA_MB_435','BT_549','T_47D'}))
	    originTissuesArray{end+1} = 'breast';
	    INITFilesArray{end+1} = 'iBreastCancer1746';
	    mCADREFilesArray{end+1} = 'breast_tumor';
	elseif any(strcmp(cellLinesArray{i},{'OVCAR_3','OVCAR_4','OVCAR_5','OVCAR_8','IGROV1','SK_OV_3'}))
	    originTissuesArray{end+1} = 'ovarian';
	    INITFilesArray{end+1} = 'iOvarianCancer1620';
	    mCADREFilesArray{end+1} = 'ovary_tumor';
	elseif any(strcmp(cellLinesArray{i},{'CCRF_CEM','K562','MOLT_4','HL_60_TB_','RPMI_8226','SR'}))
	    originTissuesArray{end+1} = 'leukemia';
	    INITFilesArray{end+1} = 'na';
	    mCADREFilesArray{end+1} = 'na';
	elseif any(strcmp(cellLinesArray{i},{'UO_31','SN12C','A498','CAKI_1','RXF_393','786_O','ACHN','TK_10'}))
	    originTissuesArray{end+1} = 'renal';
	    INITFilesArray{end+1} = 'iRenalCancer1448';
	    mCADREFilesArray{end+1} = 'kidney_tumor';
	elseif any(strcmp(cellLinesArray{i},{'LOX_IMVI','MALME_3M','SK_MEL_2','SK_MEL_5','SK_MEL_28','M14','UACC_62','UACC_257'}))
	    originTissuesArray{end+1} = 'melanoma';
	    INITFilesArray{end+1} = 'iSkinCancer1386';
	    mCADREFilesArray{end+1} = 'na';
	elseif any(strcmp(cellLinesArray{i},{'PC_3','DU_145'}))
	    originTissuesArray{end+1} = 'prostate';
	    INITFilesArray{end+1} = 'iProstateCancer1560';
	    mCADREFilesArray{end+1} = 'na';
	elseif any(strcmp(cellLinesArray{i},{'SNB_19','SNB_75','U251','SF_268','SF_295','SF_539'}))
	    originTissuesArray{end+1} = 'cns';
	    INITFilesArray{end+1} = 'na';
	    mCADREFilesArray{end+1} = 'na';
	else
	    originTissuesArray{end+1} = 'na';
	    INITFilesArray{end+1} = 'na';
	    mCADREFilesArray{end+1} = 'na';
	end
    end
end