
function cellData = predictCellTimecourse(CellLine, funToFit, otherCellLines)


%% with left-out cell lines, optimize prediction parameters
fitSettings = optimizeFittingParameters(otherCellLines, funToFit);
%fitSettings = [8 2 .1];

%% Load all data associated with current cell line
cellData = loadCellLineData(CellLine);
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
[fitParms_kp,~, ~, ~,solutions] = run(ms,prob,custpts);

% run at the solution to calculate the residual and jacobian to calculate
% the confidence intervals for each parameter
[fitParms_kp,~,res,~,~,~,J] = ...
    lsqnonlin(@(x)runGrowthModel_allData(x,cellData),fitParms_kp,...
    [], [], lsqOpts);%[0 1], [10 50000]
fitParms_kp_ci = nlparci(fitParms_kp,res,'jacobian',J);
kp = fitParms_kp(1); %hr^-1
theta = fitParms_kp(2); %cell count

%% leave one out analysis

% partition data to do leave-one-out analysis
c_part = cvpartition(cellData.ExposureTimes,'LeaveOut');

numBootstrap = 500;
%numBootstrap = 1;

%funToFit = [1 5];
numFunctionsFit = length(funToFit);

numConds_testSet = length(cellData.allYall{1})-1;

predictedTimecourses = cell(numConds_testSet, numBootstrap,c_part.NumTestSets);
predictedTimecourses_smooth = cell(numConds_testSet, numBootstrap,c_part.NumTestSets);

meanPercentError = zeros(numConds_testSet, numBootstrap, c_part.NumTestSets);
maxPercentError = zeros(numConds_testSet, numBootstrap, c_part.NumTestSets);
eoePercentError = zeros(numConds_testSet, numBootstrap, c_part.NumTestSets);

fitParms_bootstrap = cell(numFunctionsFit,numBootstrap,c_part.NumTestSets);
predictedParms_bootstrap = cell(numFunctionsFit,numBootstrap,c_part.NumTestSets);

for datasetIter = 1:c_part.NumTestSets
    
    % collect Training Data
    ioi_train = c_part.training(datasetIter);
    
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
    
    predictorVariables_training = [cb_training auc_training];
    
    %sort training data (for regularization calculation)
    [predictorVariables_training, srtedInds] = sortrows(predictorVariables_training);
    allTspan_training = allTspan_training(srtedInds);
    allTrainingData = allTrainingData(srtedInds);
    
    % collect Test Data
    ioi_test = c_part.test(datasetIter);
    
    % remove control data (no drug added)
    allTestData = cat(2,cellData.allYall{ioi_test});
    allTestData = allTestData(2:end,:);
    allTestData = allTestData(:);
    
    allTspan_test = cellData.Tspan(ioi_test);
    allTspan_test = repmat(allTspan_test,size(allTestData,1),1);
    allTspan_test = allTspan_test(:);
    
    cb_test = cat(2,cellData.dox_pss{ioi_test});
    cb_test = cb_test(2:end,:);
    cb_test = cb_test(:);
    
    auc_test = cat(2,cellData.auc{ioi_test});
    auc_test = auc_test(2:end,:);
    auc_test = auc_test(:);
    
    predictorVariables_test = [cb_test auc_test];
    
    if datasetIter == 3
        fprintf('%1d\n',cellData.ExposureTimes(ioi_test))
    else
        fprintf('%1d ',cellData.ExposureTimes(ioi_test))
    end
    
    parfor bootstrapIter = 1:numBootstrap
    %for bootstrapIter = 1:numBootstrap
        %disp(bootstrapIter)
        
        % bootstrapping to re-sample six datasets from each treatment
        % condition
        resampd_TrainingData = resampleTrainingData(allTrainingData);
        
        % within each resampled set, get fits for each treatment
        % response model
        collectAICs = zeros(size(predictorVariables_training,1),numFunctionsFit);
        hold_fitParms = cell(1, numFunctionsFit);
        
        numTestConds = size(predictorVariables_test, 1);
        hold_predictedParms = cell(1, numFunctionsFit);
        hold_modelPredictions = cell(numTestConds, numFunctionsFit);
        hold_modelPredictions_smooth = cell(numTestConds, numFunctionsFit);
        
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
                fitSettings(functIter), log(predictorVariables_training));
            
            collectAICs(:,functIter) = allAICs;
            
            hold_fitParms{functIter} = allParms;
            
            % Local regression to predict test set model parameters
            % from training model fits
            predictedParms = zeros(size(predictorVariables_test,1),numParms);
            validSpan = fitSettings(end);% try to validate this on other cell lines
            %keyboard
            for parmIter = 1:numParms
                predictedParms(:, parmIter) = ...
                    mylowess([log(predictorVariables_training), allParms(:, parmIter)],...
                    log(predictorVariables_test), validSpan);
            end
            
            if cellData.ExposureTimes(ioi_test)==12 && functIter==1 && bootstrapIter == 1 && false
                makeFigure6_localregression;
            end
            
            hold_predictedParms{functIter} = predictedParms;
            
            % Run model forward
            [allPredictions, allPredictions_smooth, smooth_timeVec] = ...
                evalForwardModel(predictedParms,...
                fixedParms, allTestData, allTspan_test, postTreat_TP);
            
            
            hold_modelPredictions(:,functIter) = allPredictions;
            hold_modelPredictions_smooth(:,functIter) = allPredictions_smooth;
            
            
        end
        
        fitParms_bootstrap(:, bootstrapIter, ioi_test) = hold_fitParms;
        predictedParms_bootstrap(:, bootstrapIter, ioi_test) = hold_predictedParms;
        %keyboard
        % fit logit function for model selection/blending
        trainVar = [log(predictorVariables_training) hold_fitParms{1,:}];
        predVar = [log(predictorVariables_test) hold_predictedParms{1,:}];
        
        modelWeights = modelWeightCalculation(collectAICs, trainVar,...
            predVar,'continuous');
        
        if cellData.ExposureTimes(ioi_test)==12 && bootstrapIter == 1 && false
            makeFigure6_weighting;
        end
        
        % blend models to make final prediction
        [blended, blended_smooth] = blendModels(modelWeights, ...
            hold_modelPredictions, hold_modelPredictions_smooth);
        
        predictedTimecourses(:,bootstrapIter,ioi_test) = blended;
        predictedTimecourses_smooth(:,bootstrapIter,ioi_test) = blended_smooth;
        
        
        % calculated normalized error
        normError = cellfun(@(x,y) abs(x-y(postTreat_TP:end,:))./y(postTreat_TP:end,:),...
            blended, allTestData,'UniformOutput',false);
        
        % average percent error over course of experiment (exclude first
        % point, since that is always correct)
        meanPE = cellfun(@(x) mean(mean(x(2:end,:))), normError);
        meanPercentError(:,bootstrapIter,ioi_test) = meanPE;
        
        % maximum percent error over course of experiment
        maxPE = cellfun(@(x) max(x(:)), normError);
        maxPercentError(:,bootstrapIter,ioi_test) = maxPE;
        
        % average percent error over final three timepoints in experiment
        meanEoE = cellfun(@(x) mean(mean(x(end-2:end,:))), normError);
        eoePercentError(:,bootstrapIter,ioi_test) = meanEoE;
        
    end
    
end

% save all results for each cell line
sfn = sprintf('%s_Predictions.mat',CellLine);
save(sfn)
