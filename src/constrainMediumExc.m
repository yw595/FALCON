function modelConstrained = constrainMediumExc(model,useMinUptake,option,excludeCOREConstrain)
% This function takes an original model, usually recon 2, and 
% constrains all models to only have very small uptake values unless
% they are one of the key medium components.

% Yiping Wang      08/??/2013
% Brandon Barker   09/21/2013    Changed to allow small amounts of
%                                uptake not listed for the medium.

if ~exist('useMinUptake','var')
    useMinUptake = 1;
end
if ~exist('option','var')
    option = 'regular';
    excludeCOREConstrain = {};
end
minUptake = -0.005;
% The above value (-0.005) reduces growth rate, but it is still far more
% than is likely in any cancer cell line.

mediumExcIdxs = loadMediumExcIdxs(model);
modelConstrained = model;
modelConstrained.description = [modelConstrained.description '_med'];


for i = 1:length(modelConstrained.rxns)
    if (~isempty(regexp(modelConstrained.rxns{i}, '^EX(.)*\(e\)$')))
        if (sum(mediumExcIdxs == i) == 0)
            if useMinUptake
                if modelConstrained.lb(i) < minUptake
                    modelConstrained.lb(i) = minUptake;
                end
            else
                modelConstrained.lb(i) = 0;
            end
        end
    end
end

if strcmp(option,'constrainCORE')
    [cellLinesArray metsArray coreTable FVAVminArray FVAVmaxArray] = readJainTable();
    coreCol = coreTable(:,strcmp(cellLinesArray,'UACC_257'));
    jainMetsToExcIdxs = loadJainMetsToExcIdxs(metsArray,model,1);
    for i=1:length(metsArray)
        if ~any(strcmp(metsArray{i},excludeCOREConstrain))
            %disp('HERE')
            excIdxs = jainMetsToExcIdxs(metsArray{i});
            if coreCol(i)<0
                for j=1:length(excIdxs)
                    modelConstrained.ub(excIdxs(j))=0;
                end
            elseif coreCol(i)>0
                for j=1:length(excIdxs)
                    modelConstrained.lb(excIdxs(j))=0;
                end
            else
                for j=1:length(excIdxs)
                    modelConstrained.ub(excIdxs(j))=0;
                    modelConstrained.lb(excIdxs(j))=0;
                end
            end
        end
    end
end

