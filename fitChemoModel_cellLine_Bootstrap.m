
function [cellData, fitParms_bootstrap, ...
    fitParms_kp, fitParms_kp_ci,...
    bestFitModels, bestFitModels_smooth] = ...
    fitChemoModel_cellLine_Bootstrap(CellLine, funToFit)

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

numFunctionsFit = length(funToFit);
numExpTime = length(cellData.ExposureTimes);

numBootstrap = 500;

numTreatConds = length(cellData.allYall{1})-1;

bestFitModels = cell(numTreatConds, numBootstrap, numExpTime);
bestFitModels_smooth = cell(numTreatConds, numBootstrap, numExpTime);

meanPercentError = zeros(numTreatConds, numBootstrap, numExpTime);
maxPercentError = zeros(numTreatConds, numBootstrap, numExpTime);
eoePercentError = zeros(numTreatConds, numBootstrap, numExpTime);

fitParms_bootstrap = cell(numFunctionsFit, numBootstrap, numExpTime);

postTreat_TP = cellData.DrugAdded_Tp + 1;

for expTimeIter = 1:numExpTime
    
    fprintf('Fits for %s, %d hr dataset\n',CellLine, cellData.ExposureTimes(expTimeIter))
    
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
    
    %fprintf([repmat('.',1,numBootstrap) '\n']);
    
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
            fixedParms = [kp theta funToFit(functIter)];
            
            % fit all concentrations at given exposure time collectively
            % no regularization
            [allParms, ~, ~, modelAIC, ~] = ...
                fitEachTreatmentModel_regularized(...
                resampd_currentData, fixedParms, postTreat_TP, allTspan,...
                0,[],resampd_weights);
            
            collectAICs(:,functIter) = modelAIC;
            hold_fitParms{functIter,1} = allParms;
            
            % run forward model evaluations for all data (not resampled
            % data)
            [eachModelFit, eachModelFit_smooth, ~] = ...
                evalForwardModel(allParms,...
                fixedParms, currentData, allTspan, postTreat_TP);
            
            hold_eachModelFit(:,functIter) = eachModelFit;
            hold_eachModelFit_smooth(:,functIter) = eachModelFit_smooth;
            
        end
        
        fitParms_bootstrap(:, bootstrapIter, expTimeIter) = hold_fitParms(:,1);
        
        minVal = min(collectAICs,[],2);
        dAIC = collectAICs - repmat(minVal, 1, size(collectAICs,2));
        expdAIC = exp(-0.5.*dAIC);
        modelWeights = expdAIC ./ repmat(sum(expdAIC,2),1,size(expdAIC,2));
        
        [blended, blended_smooth] = blendModels(modelWeights, ...
            hold_eachModelFit, hold_eachModelFit_smooth);
        
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

save(fullfile(pwd, ['BootstrapFits_IJcounts_' CellLine '.mat']))
