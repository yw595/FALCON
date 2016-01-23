function nci60ScriptHelper2(outputFile1, outputFile2, v_fbaIn, modelToRun)
    v_fba = v_fbaIn;
    save(outputFile1, 'v_fba');
    fluxOutputFI = fopen(outputFile2,'w');
    for l=1:length(v_fba)
        fprintf(fluxOutputFI,['R_' modelToRun.rxns{l} ',' num2str(v_fba(l)) '\n']);
    end
    fclose(fluxOutputFI);
end