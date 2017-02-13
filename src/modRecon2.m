function newRecon2 = modRecon2(origRecon2)

newRecon2 = origRecon2;
newRecon2.description = 'newRecon2';
for i=1:length(newRecon2.grRules)
    disp(['processing gene ' num2str(i)])
    %check if gene name contains decimal
    if ~isempty(regexp(newRecon2.grRules{i},'\.'))
        x = newRecon2.grRules{i};
        %strip out decimal and digits after decimal,
        %stopping after EITHER a right parenthesis or space is encountered
        while ~isempty(regexp(x,'\.\d+(\)| )'))
            startIdxs = regexp(x,'\.\d+(\)| )','start');
            endIdxs = regexp(x,'\.\d+(\)| )','end');
            x = [x(1:startIdxs(1)-1) x(endIdxs(1):end)];
        end
        %strip out decimal and digits after decimal if they run to end of gene id
        if ~isempty(regexp(x,'\.\d+$'))
            x = x(1:regexp(x,'\.')-1);
        end
        newRecon2.grRules{i} = x;
    end
end