
function MPE = calculate_MeanPercentError(...
    trainParms, collectAICs, predictorVariables_training,...
    cell_kpTheta, funToFit, parmSpanVals,...
    allTestData, allTspan_test, predictorVariables_test,...
    postTreat_TP,weightCell)


hold_predictedParms = cell(1, length(funToFit));
hold_modelPredictions = cell(size(allTestData,1), length(funToFit));

for functIter = 1:length(funToFit)
    
    allParms = trainParms{functIter};
    fixedParms = [cell_kpTheta funToFit(functIter)];
    
    % Local regression to predict test set model parameters
    % from training model fits
    numParms = size(allParms,2);
    predictedParms = zeros(length(predictorVariables_test),numParms);
    
    parfor parmIter = 1:numParms
    %for parmIter = 1:numParms
        
        c_validSpan = parmSpanVals;
        
        predictedParms(:, parmIter) = ...
            mylowess([log(predictorVariables_training), allParms(:, parmIter)],...
            log(predictorVariables_test), c_validSpan);
        
    end
    
    
    hold_predictedParms{functIter} = predictedParms;
    
    % Run model forward
    allPredictions = ...
        evalForwardModel(predictedParms,...
        fixedParms, allTestData, allTspan_test, postTreat_TP);
    
    hold_modelPredictions(:,functIter) = allPredictions;
    
end

% fit logit function for model selection/blending
trainVar = [log(predictorVariables_training) trainParms{1,:}];
predVar = [log(predictorVariables_test) hold_predictedParms{1,:}];

modelWeights = modelWeightCalculation(collectAICs, trainVar,...
    predVar,'continuous');

% blend models to make final prediction
blended = blendModels(modelWeights, ...
    hold_modelPredictions, hold_modelPredictions);

% calculated normalized error
normError = cellfun(@(x,y,z) ...
    (abs(x-y(postTreat_TP:end,:))./y(postTreat_TP:end,:)) .* z(postTreat_TP:end,:),...
    blended, allTestData, weightCell,'UniformOutput',false);

% average percent error over course of experiment (exclude first
% point, since that is always correct)
meanPE_cell = cellfun(@(x) x(2:end,:),normError,'UniformOutput',false);
MPE = mean(cellfun(@(x) mean(mean(x(x>0))),meanPE_cell));


