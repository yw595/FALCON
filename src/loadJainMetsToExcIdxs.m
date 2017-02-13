function jainMetsToExcIdxs = loadJainMetsToExcIdxs(jainMetsArray,model,useReducedNames)
if ~exist('useReducedNames','var')
    useReducedNames = 0;
end
jainMetsToExcIdxs=containers.Map;
for i=1:length(jainMetsArray)
    if(strcmp(jainMetsArray{i},'34hpp'))
        excrxnname='EX_34hpp';
        excrxnind=find(strcmp(excrxnname,model.rxns));
        jainMetsToExcIdxs(jainMetsArray{i})=excrxnind;
    elseif(strcmp(jainMetsArray{i},'glc_D'))
        excrxnname1='EX_glc(e)';
        excrxnind1=find(strcmp(excrxnname1,model.rxns));
        if useReducedNames
            excrxnname2='TR_glc_D[c]';
            excrxnind2=find(strcmp(excrxnname2,model.rxns));
            jainMetsToExcIdxs(jainMetsArray{i})=[excrxnind1 excrxnind2];
        else
            jainMetsToExcIdxs(jainMetsArray{i})=excrxnind1;
        end
    elseif(strcmp(jainMetsArray{i},'udpgal/udpg'))
        excrxnname1='EX_udpgal(e)';
        excrxnind1=find(strcmp(excrxnname1,model.rxns));
        excrxnname2='EX_udpg(e)';
        excrxnind2=find(strcmp(excrxnname2,model.rxns));
        jainMetsToExcIdxs(jainMetsArray{i})=[excrxnind1 excrxnind2];
    elseif(strcmp(jainMetsArray{i},'lac_D/lac_L'))
        excrxnname1='EX_lac_D(e)';
        excrxnind1=find(strcmp(excrxnname1,model.rxns));
        excrxnname2='EX_lac_L(e)';
        excrxnind2=find(strcmp(excrxnname2,model.rxns));
        if useReducedNames
            excrxnname3='TR_lac_L[c]';
            excrxnind3=find(strcmp(excrxnname3,model.rxns));
            jainMetsToExcIdxs(jainMetsArray{i})=[excrxnind1 excrxnind2 excrxnind3];
        else
            jainMetsToExcIdxs(jainMetsArray{i})=[excrxnind1 excrxnind2];
        end
    elseif(strcmp(jainMetsArray{i},'sbt_D'))
        excrxnname='EX_sbt-d(e)';
        excrxnind=find(strcmp(excrxnname,model.rxns));
        jainMetsToExcIdxs(jainMetsArray{i})=excrxnind;
    elseif(strcmp(jainMetsArray{i},'tyr_l'))
        excrxnname1='EX_tyr_L(e)';
        excrxnind1=find(strcmp(excrxnname1,model.rxns));
        if useReducedNames
            excrxnname2='TR_tyr_L[c]';
            excrxnind2=find(strcmp(excrxnname2,model.rxns));
            excrxnname3='DEMAND_tyr_L[c]';
            excrxnind3=find(strcmp(excrxnname3,model.rxns));
            jainMetsToExcIdxs(jainMetsArray{i})=[excrxnind1 excrxnind2 excrxnind3];
        else
            jainMetsToExcIdxs(jainMetsArray{i})=excrxnind1;
        end
    elseif(strcmp(jainMetsArray{i},'cit/icit'))
        excrxnname='EX_cit(e)';
        excrxnind=find(strcmp(excrxnname,model.rxns));
        jainMetsToExcIdxs(jainMetsArray{i})=excrxnind;
    elseif(strcmp(jainMetsArray{i},'N/A'))
        excrxnname='';
        excrxnind=0;
        jainMetsToExcIdxs(jainMetsArray{i})=excrxnind;
    else
        excrxnname1=strcat('EX_',strcat(jainMetsArray{i},'(e)'));
        excrxnind1=find(strcmp(excrxnname1,model.rxns));
        if useReducedNames
            excrxnindarr = [excrxnind1];
            excrxnname2=strcat('TR_',strcat(jainMetsArray{i},'[c]'));
            excrxnind2=find(strcmp(excrxnname2,model.rxns));
            if ~isempty(excrxnind2)
                excrxnindarr(end+1) = excrxnind2;
            end
            excrxnname3=strcat('DEMAND_',strcat(jainMetsArray{i},'[c]'));
            excrxnind3=find(strcmp(excrxnname3,model.rxns));
            if ~isempty(excrxnind3)
                excrxnindarr(end+1) = excrxnind3;
            end
            if length(excrxnindarr)==1
                excrxnindarr = excrxnindarr(1);
            end
            jainMetsToExcIdxs(jainMetsArray{i})=excrxnindarr;
        else
            jainMetsToExcIdxs(jainMetsArray{i})=excrxnind1;
        end
    end
end

if useReducedNames
end
end

