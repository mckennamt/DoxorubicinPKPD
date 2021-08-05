
function cellData = predictCellTimecourse_only12(CellLine, funToFit, otherCellLines)

% function used to train model with only 12-hour exposure data to predict
% response to both 6-hour and 24-hour exposures

%% with left-out cell lines, optimize prediction parameters
%fitSettings = optimizeFittingParameters(otherCellLines, funToFit);
fitSettings = optimizeFittingParameters_only12(otherCellLines, funToFit);
%fitSettings = [5 3 .2];

%% Load all data associated with current cell line
%cellData = loadCellLineData(CellLine);
cellData = loadCellLineData_IJ(CellLine);
cellData = rmfield(cellData,{'drugConc','NumNuclei','NumNuclei2',...
    'NucDensity','DrugTime','times'});
fprintf('%s\n',CellLine)

%% Get fits for growth rate, theta
% Fit to control data and pre-treatment timepoints
% might move to prediction section

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
fitParms_kp = run(ms,prob,custpts);

% run at the solution to calculate the residual and jacobian to calculate
% the confidence intervals for each parameter
[fitParms_kp,~,res,~,~,~,J] = ...
    lsqnonlin(@(x)runGrowthModel_allData(x,cellData),fitParms_kp,...
    [], [], lsqOpts);%[0 1], [10 50000]
fitParms_kp_ci = nlparci(fitParms_kp,res,'jacobian',J);
kp = fitParms_kp(1); %hr^-1
theta = fitParms_kp(2); %cell count

%% train with 12 hour data to predict 6 and 24 hour data sets
% use only 12 hour data to predict 6 hour and 24 hour sets

% bootstrap sample training data to generate CI for predictions
numBootstrap = 500;
%numBootstrap = 1;

%funToFit = [1 5];
numFunctionsFit = length(funToFit);

numConds_testSet = length(cellData.allYall{1})-1; %includes a control

predictedTimecourses = cell(numConds_testSet, numBootstrap,2);
predictedTimecourses_smooth = cell(numConds_testSet, numBootstrap,2);

meanPercentError = zeros(numConds_testSet, numBootstrap, 2);
maxPercentError = zeros(numConds_testSet, numBootstrap, 2);
eoePercentError = zeros(numConds_testSet, numBootstrap, 2);

fitParms_bootstrap = cell(numFunctionsFit,numBootstrap,1);% only fitting 12 hour data
predictedParms_bootstrap = cell(numFunctionsFit,numBootstrap,2);


% collect Training Data
ioi_train = logical([0 1 0]');

% remove control data (no drug added)
allTrainingData = cat(2,cellData.allYall{ioi_train});
allTrainingData = allTrainingData(2:end,:);
numConcs = size(allTrainingData,1);
allTrainingData = allTrainingData(:);

allTspan_training = cellData.Tspan(ioi_train);
allTspan_training = repmat(allTspan_training,numConcs,1);
allTspan_training = allTspan_training(:);

cb_training = cat(2,cellData.dox_pss{ioi_train});
cb_training = cb_training(2:end,:);
cb_training = cb_training(:);

auc_training = cat(2,cellData.auc{ioi_train});
auc_training = auc_training(2:end,:);
auc_training = auc_training(:);

% extract weights
if isfield(cellData, 'allYweights')
    weights_training = cat(2,cellData.allYweights{ioi_train});
    weights_training = weights_training(2:end,:);
    weights_training = weights_training(:);
else
    weights_training = [];
end

predictorVariables_training = [cb_training auc_training];

%sort training data (for regularization calculation)
[predictorVariables_training, srtedInds] = sortrows(predictorVariables_training);
allTspan_training = allTspan_training(srtedInds);
allTrainingData = allTrainingData(srtedInds);
weights_training = weights_training(srtedInds);

% collect Test Data (6 hour, 24 hour)
ioi_testGroups = logical([1 0 0; 0 0 1]);
allTestData = cell(1,size(ioi_testGroups,1));
allTspan_test = cell(1,size(ioi_testGroups,1));
predictorVariables_test = cell(1,size(ioi_testGroups,1));
allTestWeights = cell(1,size(ioi_testGroups,1));

for tgi = 1:size(ioi_testGroups,1)
    
    ioi_test = ioi_testGroups(tgi,:);
    
    % remove control data (no drug added)
    allTestData_tg = cat(2,cellData.allYall{ioi_test});
    allTestData_tg = allTestData_tg(2:end,:);
    allTestData{tgi} = allTestData_tg(:);
    
    allTspan_test_tg = cellData.Tspan(ioi_test);
    allTspan_test_tg = repmat(allTspan_test_tg,size(allTestData_tg,1),1);
    allTspan_test{tgi} = allTspan_test_tg(:);
    
    cb_test = cat(2,cellData.dox_pss{ioi_test});
    cb_test = cb_test(2:end,:);
    cb_test = cb_test(:);
    
    auc_test = cat(2,cellData.auc{ioi_test});
    auc_test = auc_test(2:end,:);
    auc_test = auc_test(:);
    
    if isfield(cellData, 'allYweights')
        weights_test = cat(2,cellData.allYweights{ioi_test});
        weights_test = weights_test(2:end,:);
        weights_test = weights_test(:);
    else
        weights_test = [];
    end
    allTestWeights{tgi} = weights_test;
    
    predictorVariables_test{tgi} = [cb_test auc_test];
end

%parfor bootstrapIter = 1:numBootstrap
for bootstrapIter = 1:numBootstrap
    
    fprintf('%2d\n',bootstrapIter)
    
    % bootstrapping to re-sample six datasets from each treatment
    % condition
    [resampd_TrainingData,resampd_weights] = ...
        resampleTrainingData(allTrainingData,weights_training);
    
    % within each resampled set, get fits for each treatment
    % response model
    collectAICs = zeros(size(predictorVariables_training,1),numFunctionsFit);
    hold_fitParms = cell(1,numFunctionsFit);
    
    numTestConds = size(predictorVariables_test{1}, 1);
    hold_predictedParms = cell(numFunctionsFit,1,2);
    hold_modelPredictions = cell(numTestConds, numFunctionsFit,2);
    hold_modelPredictions_smooth = cell(numTestConds, numFunctionsFit,2);
    
    for functIter = 1:numFunctionsFit
        
        % fit current model
        fixedParms = [kp theta funToFit(functIter)];
        postTreat_TP = cellData.DrugAdded_Tp + 1;
        
        %[allParms, allParms_ci, allInitGuess, allAICs, numParms] = ...
        %    fitEachTreatmentModel(...
        %    resampd_TrainingData, fixedParms, postTreat_TP, allTspan_training);
        [allParms, allParms_ci, allInitGuess, allAICs, numParms] = ...
            fitEachTreatmentModel_regularized(...
            resampd_TrainingData, fixedParms, postTreat_TP, allTspan_training,...
            fitSettings(functIter), log(predictorVariables_training),...
            resampd_weights);
        
        collectAICs(:,functIter) = allAICs;
        hold_fitParms{functIter} = allParms;
        
        % generate parameter estimates for each test set
        % Local regression to predict test set model parameters
        % from training model fits
        for testCond = 1:size(ioi_testGroups,1)
            predictedParms = zeros(size(predictorVariables_test{testCond},1),numParms);
            validSpan = fitSettings(end);% try to validate this on other cell lines
            for parmIter = 1:numParms
                predictedParms(:, parmIter) = ...
                    mylowess([log(predictorVariables_training), allParms(:, parmIter)],...
                    log(predictorVariables_test{testCond}), validSpan);
            end
            hold_predictedParms{functIter,1,testCond} = predictedParms;
            
            % Run model forward
            [allPredictions, allPredictions_smooth, smooth_timeVec] = ...
                evalForwardModel(predictedParms,...
                fixedParms, allTestData{testCond},...
                allTspan_test{testCond}, postTreat_TP);
            hold_modelPredictions(:,functIter,testCond) = allPredictions;
            hold_modelPredictions_smooth(:,functIter,testCond) = allPredictions_smooth;
            
        end
        
        if cellData.ExposureTimes(2)==12 && functIter==1 && bootstrapIter == 1 && false
            makeFigure6_localregression;
        end
        
    end
    
    fitParms_bootstrap(:, bootstrapIter) = hold_fitParms;
    predictedParms_bootstrap(:, bootstrapIter, :) = hold_predictedParms;
    
    if cellData.ExposureTimes(2)==12 && bootstrapIter == 1 && false
        makeFigure6_weighting;
    end
    
    % fit logit function for model selection/blending
    trainVar = [log(predictorVariables_training) hold_fitParms{1,:}];
    
    % define these to allow for parallel running
    hold_predictedTimecourses = cell(numTestConds,2);
    hold_predictedTimecourses_smooth = cell(numTestConds,2);
    hold_meanPercentError = zeros(numTestConds,2);
    hold_maxPercentError = zeros(numTestConds,2);
    hold_eoePercentError = zeros(numTestConds,2);
    
    for testCond = 1:size(ioi_testGroups,1)
        pp_training = hold_predictedParms(:,1,testCond)';
        predVar = [log(predictorVariables_test{testCond}) pp_training{1,:}];
        
        modelWeights = modelWeightCalculation(collectAICs, trainVar,...
            predVar,'continuous');
        
        % blend models to make final prediction, 6 hours
        [blended, blended_smooth] = blendModels(modelWeights, ...
            hold_modelPredictions(:,:,testCond),...
            hold_modelPredictions_smooth(:,:,testCond));
        
        %predictedTimecourses(:,bootstrapIter,testCond) = blended;
        %predictedTimecourses_smooth(:,bootstrapIter,testCond) = blended_smooth;
        hold_predictedTimecourses(:,testCond) = blended;
        hold_predictedTimecourses_smooth(:,testCond) = blended_smooth;
        
        % calculated normalized error
        cWeights = allTestWeights{testCond};
        if isempty(cWeights)
            weightCell = cellfun(@(x) ones(size(x)), allTestData{testCond}, 'UniformOutput',false);
        else
            weightCell = cellfun(@(x) zeros(size(x)), allTestData{testCond}, 'UniformOutput',false);
            for cellIter = 1:length(weightCell)
                for colIter = 1:size(weightCell{cellIter},2)
                    weightCell{cellIter}(1:cWeights{cellIter}(colIter),colIter) = 1;
                end
            end
        end
        normError = cellfun(@(x,y,z) ...
            (abs(x-y(postTreat_TP:end,:))./y(postTreat_TP:end,:)) .* z(postTreat_TP:end,:),...
            blended, allTestData{testCond}, weightCell,'UniformOutput',false);
        
        % average percent error over course of experiment (exclude first
        % point, since that is always correct)
        meanPE_nz = cellfun(@(x,y) x(2:end,:),...
            normError,'UniformOutput',false);
        meanPE = cellfun(@(x) mean(mean(x(x>0))),meanPE_nz);
        %meanPercentError(:,bootstrapIter,testCond) = meanPE;
        hold_meanPercentError(:,testCond) = meanPE;
        
        % maximum percent error over course of experiment
        maxPE = cellfun(@(x) max(x(:)), normError);
        %maxPercentError(:,bootstrapIter,testCond) = maxPE;
        hold_maxPercentError(:,testCond) = maxPE;
        
        % average percent error over final three timepoints in experiment
        meanEoE = zeros(length(meanPE_nz),1);
        for cellIter = 1:length(meanPE_nz)
            colErrors = zeros(1,size(meanPE_nz{cellIter},2));
            for colIter = 1:size(meanPE_nz{cellIter},2)
                hold_nz = meanPE_nz{cellIter}(logical(weightCell{cellIter}(postTreat_TP+1:end,colIter)),colIter);
                colErrors(colIter) = mean(hold_nz(end-2:end));
            end
            meanEoE(cellIter) = mean(colErrors);
        end
        %meanEoE = cellfun(@(x) mean(mean(x(end-2:end,:))), normError);
        %eoePercentError(:,bootstrapIter,testCond) = meanEoE;
        hold_eoePercentError(:,testCond) = meanEoE;
        
    end
    predictedTimecourses(:,bootstrapIter,:) = hold_predictedTimecourses;
    predictedTimecourses_smooth(:,bootstrapIter,:) = hold_predictedTimecourses_smooth;
    meanPercentError(:,bootstrapIter,:) = hold_meanPercentError;
    maxPercentError(:,bootstrapIter,:) = hold_maxPercentError;
    eoePercentError(:,bootstrapIter,:) = hold_eoePercentError;
    
    
end


% save all results for each cell line
sfn = sprintf('%s_Predictions_IJCounts_only12.mat',CellLine);
save(sfn)
