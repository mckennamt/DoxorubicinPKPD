% get PD parameter values from all data
% main function to call all fitting fuctions

rng('default');
addpath('./Data');

% Create parallel pool
poolobj = gcp('nocreate');
if isempty(poolobj)
    c = parcluster('local');
    poolsize = c.NumWorkers;
    myPool = parpool(poolsize);
else
    myPool = poolobj;
end

% Specify data to be used
allCellLine = {'SUM149','MDAMB468', 'MDAMB231', 'MDAMB453'};
numberCellLines = length(allCellLine);
% Published functions (apoptosis (1), mitotic catastrophe (5))
funToFit = [1 5];

allCellData = cell(size(allCellLine));
collect_Fits = cell(size(allCellLine));
collect_kp = cell(size(allCellLine));
collect_AIC = cell(size(allCellLine));
collect_FitsCI = cell(size(allCellLine));

for cellLineIter = 1:length(allCellLine)
    
    CellLine = allCellLine{cellLineIter};
    
    %[allCellData{cellLineIter}, collect_AIC{cellLineIter},...
    %    collect_Fits{cellLineIter}, collect_FitsCI{cellLineIter},...
    %    collect_kp{cellLineIter}, ~] = ...
    %    fitChemoModel_cellLine(CellLine, funToFit, [1 5]); 
    
    %[~,~,~,~,~,~] = ...
    %    fitChemoModel_cellLine_Bootstrap(CellLine, funToFit); 
    
    [~,~,~,~,~,~] = ...
        fitGrowthModel_cellLine_Bootstrap(CellLine);
    
    
end

delete(myPool)