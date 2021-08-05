
function [cellData, collectAICs, parameterFits, ...
    parameterCI, fitParms_kp, fitParms_kp_ci] = ...
    fitChemoModel_cellLine(CellLine, funToFit, regConstant)

% function to fit treatment response data to each cell line
% this function does not use a bootstrap approach to fit data, simply fits
% the six replicates for each treatment condition simultaneously

%% load data corresponding to specified cell line
%cellData = loadCellLineData(CellLine);
cellData = loadCellLineData_IJ(CellLine);
cellData = rmfield(cellData,{'drugConc','NumNuclei','NumNuclei2',...
    'NucDensity','DrugTime','times'});

%% Get fits for growth rate, theta (carrying capacity)
% Fit to control data and pre-treatment timepoints
% use all exposure time data sets for this
lsqOpts = optimoptions('lsqnonlin','MaxFunEvals',10000,'Display','none',...
    'TolFun', 1e-50,'MaxIter', 5000,'TolX', 1e-5);
prob = createOptimProblem('lsqnonlin','x0',[0 0],...
    'objective',@(x) runGrowthModel_allData(x,cellData),...
    'lb',[0 1],'ub',[1 1e6],'options',lsqOpts);

% initialize with 25 sets of parameter estimates to avoid local minima
ms = MultiStart('UseParallel',true,'Display','off');

% create parameter estimates
parmEsts{1} = linspace(.01, .05, 5);
parmEsts{2} = linspace(5000, 1e5, 5);
hldGrid = cell(1,numel(parmEsts));
[hldGrid{:}] = ndgrid(parmEsts{:});
hldGrid = cellfun(@(x) x(:), hldGrid,'UniformOutput',false);
startPointList = cat(2,hldGrid{:});
custpts = CustomStartPointSet(startPointList);
[fitParms_kp,~, ~, ~,ms_solutions] = run(ms,prob,custpts);

% run at the solution to calculate the residual and jacobian to calculate
% the confidence intervals for each parameter
[fitParms_kp,~,res,~,~,~,J] = ...
    lsqnonlin(@(x)runGrowthModel_allData(x,cellData),fitParms_kp,...
    [], [], lsqOpts);
fitParms_kp_ci = nlparci(fitParms_kp,res,'jacobian',J);
kp = fitParms_kp(1); %hr^-1
theta = fitParms_kp(2); %cell count

%% fit treatment response model

numFunctionsFit = length(funToFit);
numExpTime = length(cellData.ExposureTimes);

collectAICs = cell(numExpTime,numFunctionsFit);
parameterFits = cell(numExpTime,numFunctionsFit);
parameterCI = cell(numExpTime,numFunctionsFit);

for expTimeIter = 1:numExpTime
    
    % remove control data (no drug added)
    currentData = cellData.allYall{expTimeIter};
    currentData = currentData(2:end,:);
    
    % extract weights
    if isfield(cellData, 'allYweights')
       currentWeights = cellData.allYweights{expTimeIter}; 
       currentWeights = currentWeights(2:end,:);
    else
        currentWeights = [];
    end
    
    allTspan = cellData.Tspan(expTimeIter);
    allTspan = repmat(allTspan,size(currentData,1),1);
    
    currDoxPss = cellData.dox_pss{expTimeIter};
    currDoxPss = currDoxPss(2:end,:);
    
    currAUC = cellData.auc{expTimeIter};
    currDoxPss = [currDoxPss currAUC(2:end,:)];
    
    for functIter = 1:numFunctionsFit
        
        % fit current model
        fixedParms = [kp theta funToFit(functIter)];
        postTreat_TP = cellData.DrugAdded_Tp + 1;
        
        %%fit each concentration individually
        %[allParms, allParms_ci, allInitGuess, allAICs, numParms] = ...
        %    fitEachTreatmentModel(...
        %    currentData, fixedParms, postTreat_TP, allTspan);
        
        % fit all concentrations at given exposure time collectively
        %[allParms, allParms_ci, allInitGuess, allAICs, numParms] = ...
        %    fitEachTreatmentModel_regularized(...
        %    currentData, fixedParms, postTreat_TP, allTspan, 10)
        
        %[allParms, allParms_ci, ~, allAICs, numParms] = ...
        %    fitEachTreatmentModel_regularized(...
        %    currentData, fixedParms, postTreat_TP, allTspan, 0,[]);
        
        [allParms, allParms_ci, ~, allAICs, numParms] = ...
            fitEachTreatmentModel_regularized(...
            currentData, fixedParms, postTreat_TP, allTspan,...
            regConstant(functIter), log(currDoxPss),currentWeights);
        
        collectAICs{expTimeIter,functIter} = allAICs;
        parameterFits{expTimeIter,functIter} = allParms;
        parameterCI{expTimeIter,functIter} = allParms_ci;
        
    end
    
    
end
