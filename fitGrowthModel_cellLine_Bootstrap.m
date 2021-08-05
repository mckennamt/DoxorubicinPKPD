
function [cellData, fitParms_bootstrap, fitParms_kp, fitParms_kp_ci,...
    bestFitModels, bestFitModels_smooth] = ...
    fitGrowthModel_cellLine_Bootstrap(CellLine)

% function to fit treatment response data to each cell line
% used to generate confidence intervals for fits

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

numFunctionsFit = 1;
numExpTime = length(cellData.ExposureTimes);

numBootstrap = 4;

numTreatConds = 1;

bestFitModels = cell(numTreatConds, numBootstrap, numExpTime);
bestFitModels_smooth = cell(numTreatConds, numBootstrap, numExpTime);

meanPercentError = zeros(numTreatConds, numBootstrap, numExpTime);
maxPercentError = zeros(numTreatConds, numBootstrap, numExpTime);
eoePercentError = zeros(numTreatConds, numBootstrap, numExpTime);

fitParms_bootstrap = cell(numFunctionsFit, numBootstrap, numExpTime);

postTreat_TP = cellData.DrugAdded_Tp + 1;

for expTimeIter = 1:numExpTime
    
    fprintf('Fits for %s, %d hr dataset\n',CellLine, cellData.ExposureTimes(expTimeIter))
    
    % isolate control data (no drug added)
    currentData = cellData.allYall{expTimeIter};
    currentData = currentData(1,:);
    
    % extract weights
    if isfield(cellData, 'allYweights')
        currentWeights = cellData.allYweights{expTimeIter};
        currentWeights = currentWeights(1,:);
    else
        currentWeights = [];
    end
    
    allTspan = cellData.Tspan(expTimeIter);
    allTspan = repmat(allTspan,size(currentData,1),1);
    
    % run through fitting process 500 times to generate confidence intervals
    for bootstrapIter = 1:numBootstrap
        
        fprintf('%d ',bootstrapIter)
        %disp(bootstrapIter)
        %fprintf('\b|\n')
        
        % resample within each concentration
        [resampd_currentData, resampd_weights] = ...
            resampleTrainingData(currentData,currentWeights);
        
        hold_eachModelFit = cell(numTreatConds, numFunctionsFit);
        hold_eachModelFit_smooth = cell(numTreatConds, numFunctionsFit);
        
        collectAICs = zeros(numTreatConds,numFunctionsFit);
        hold_fitParms = cell(numFunctionsFit,1);
        
        for functIter = 1:numFunctionsFit
            
            % fit current model
            fixedParms = [kp theta];
            
            % fit all concentrations at given exposure time collectively
            % no regularization
            lsqOpts = optimoptions('lsqnonlin','MaxFunEvals',10000,'Display','none',...
                'TolFun', 1e-50,'MaxIter', 5000,'TolX', 1e-5);
            prob = createOptimProblem('lsqnonlin','x0',[0 0],...
                'objective',@(x) runGrowthModel_controlData(x,cellData),...
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
            [allParms,~, ~, ~,ms_solutions] = run(ms,prob,custpts);
            
            
            collectAICs(:,functIter) = 0;
            hold_fitParms{functIter,1} = allParms;
            
            % run forward model evaluations for all data (not resampled
            % data)
            
            [eachModelFit, eachModelFit_smooth] = ...
                evalGrowthModel_analytic(currentData, allTspan, allParms, postTreat_TP);
            
            hold_eachModelFit(:,functIter) = eachModelFit;
            hold_eachModelFit_smooth(:,functIter) = eachModelFit_smooth;
            
        end
        
        fitParms_bootstrap(:, bootstrapIter, expTimeIter) = hold_fitParms(:,1);
        
        blended = hold_eachModelFit;
        blended_smooth = hold_eachModelFit_smooth;
        
        bestFitModels(:,bootstrapIter,expTimeIter) = blended;
        bestFitModels_smooth(:,bootstrapIter,expTimeIter) = blended_smooth;
        
        % calculated normalized error
        if isempty(currentWeights)
            weightCell = cellfun(@(x) ones(size(x)), currentData, 'UniformOutput',false);
        else
            weightCell = cellfun(@(x) zeros(size(x)), currentData, 'UniformOutput',false);
            for cellIter = 1:length(weightCell)
                for colIter = 1:size(weightCell{cellIter},2)
                    weightCell{cellIter}(1:currentWeights{cellIter}(colIter),colIter) = 1;
                end
            end
        end
        normError = cellfun(@(x,y,z) ...
            (abs(x-y(postTreat_TP:end,:))./y(postTreat_TP:end,:)) .* z(postTreat_TP:end,:),...
            blended, currentData, weightCell,'UniformOutput',false);
        
        % average percent error over course of experiment (exclude first
        % point, since that is always correct)
        meanPE_cell = cellfun(@(x) x(2:end,:),normError,'UniformOutput',false);
        meanPE = cellfun(@(x) mean(mean(x(x>0))),meanPE_cell);
        meanPercentError(:,bootstrapIter,expTimeIter) = meanPE;
        
        % maximum percent error over course of experiment
        maxPE = cellfun(@(x) max(x(:)), normError);
        maxPercentError(:,bootstrapIter,expTimeIter) = maxPE;
        
        % average percent error over final three timepoints in experiment
        % only consider valid timepoints
        meanEoE = zeros(length(meanPE_cell),1);
        for cellIter = 1:length(meanPE_cell)
            colErrors = zeros(1,size(meanPE_cell{cellIter},2));
            for colIter = 1:size(meanPE_cell{cellIter},2)
                hold_nz = meanPE_cell{cellIter}(logical(weightCell{cellIter}(postTreat_TP+1:end,colIter)),colIter);
                colErrors(colIter) = mean(hold_nz(end-2:end));
            end
            meanEoE(cellIter) = mean(colErrors);
        end
        eoePercentError(:,bootstrapIter,expTimeIter) = meanEoE;
        
    end
    fprintf('\n')
    
end

save(fullfile(pwd, ['BootstrapFits_ControlData_IJcounts_' CellLine '_20170325.mat']))

end

function [allPredictions_cell,allPredictions_smooth_cell] =...
    evalGrowthModel_analytic(currentData, allTspan, fitParms, postTreat_TP)


cData = currentData{1};
Tspan = allTspan{1};
Tspan = Tspan(postTreat_TP:end);
Tspan_mod = linspace(Tspan(1),Tspan(end),100);

allPredictions = zeros(length(Tspan), size(cData,2));
allPredictions_smooth = zeros(length(Tspan_mod), size(cData,2));

for initConditIter = 1:size(cData,2)
    
    N0 = cData(postTreat_TP,initConditIter);
    
    % same timepoints as measured data
    kp = fitParms(1);
    theta = fitParms(2);
    t = (Tspan-min(Tspan));
    N = theta*N0 ./ (N0 + (theta-N0).*exp(-kp.*t));
    allPredictions(:,initConditIter) = N;
    
    %smoothed version
    t = (Tspan_mod-min(Tspan_mod));
    N = theta*N0 ./ (N0 + (theta-N0).*exp(-kp.*t));
    allPredictions_smooth(:,initConditIter) = N;
    
end

allPredictions_cell{1} = allPredictions;
allPredictions_smooth_cell{1} = allPredictions_smooth;


end
