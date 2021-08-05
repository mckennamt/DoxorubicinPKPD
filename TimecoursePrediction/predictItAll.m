% make predictions on all cell lines and generate prediction statistics

rng('default');
addpath('./Data')

poolobj = gcp('nocreate'); % If no pool, do not create new one.
if isempty(poolobj)
    c = parcluster('local');
    poolsize = c.NumWorkers;
    myPool = parpool(poolsize);
else
    myPool = poolobj;
end

allCellLine = {'SUM149','MDAMB231','MDAMB453','MDAMB468'};
numberCellLines = length(allCellLine);
funToFit = [1 5];

allCellData = cell(size(allCellLine));

pctRunOnAll warning('off')

for cellLineIter = 1%:length(allCellLine)
    
    CellLine = allCellLine{cellLineIter};
    otherCellLines = allCellLine;
    otherCellLines(cellLineIter) = [];
    
    %allCellData{cellLineIter} = ...
    %    predictCellTimecourse(CellLine, funToFit, otherCellLines);
    
    allCellData{cellLineIter} = ...
        predictCellTimecourse_only12(CellLine, funToFit, otherCellLines);
    
end

delete(myPool)