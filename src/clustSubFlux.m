function subToAbsFlux = clustSubFlux(model,dist)

%sortAbsFlux = [];
uniqSubs = unique(model.subSystems);
for i=1:length(uniqSubs)
    uniqSubs{i} = strrep(uniqSubs{i},'/','_');
    uniqSubs{i} = strrep(uniqSubs{i},' ','_');
    uniqSubs{i} = strrep(uniqSubs{i},',','_');
end
subToAbsFlux = containers.Map;
for i=1:length(uniqSubs)
    subToAbsFlux(uniqSubs{i}) = sum(abs(dist(strcmp(model.subSystems,uniqSubs{i}))));
end
%[sortAbsFlux sortIdxs] = sort(sortAbsFlux,'descend');
%sortSubs = uniqSubs(sortIdxs);

end




    