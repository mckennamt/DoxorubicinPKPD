

function modelWeights = modelWeightCalculation(AICs, trainVar, predVar,weightingScheme)

% calculate appropriate weights for each model
% AICs are calculated Akaike Information Criterion for each model
% trainVar is [p,n] matrix of p observations and n variables to train a
% model to predict weighting vectors
% predVar is [p,n] matrix to generate weighting for predicted model
% weighting scheme is used to select technique for model weight
% calculation, 'binary' to train logistic regression model or 'continuous' to use
% AIC values to calculate weights

modelWeights = zeros(size(predVar,1), size(AICs,2));

if size(AICs,2) == 1
    modelWeights(:) = 1;
    return
    
elseif size(AICs,2) > 2
    % Only working with 2 models, best AIC wins if working with more than 2
    % functions
    
    % currently defaults to model 1 as always correct, probably could
    % implement some multi-dimensional logit, but don't have to
    modelWeights(:,1) = 1;
    
    return
    
end

modelWeights = zeros(size(predVar,1), size(AICs,2));

switch lower(weightingScheme)
    
    case 'binary'
        
        
        % true model = 0 when function in first column is favored
        % true model = 1 when function in second column is favored
        [~,trueModel] = min(AICs,[],2);
        trueModel = trueModel-1;
        
        % train logistic model with input training variables (trainVar) and true
        % model classification (trueModel)
        b_coeff = glmfit(trainVar, [trueModel ones(size(trueModel))],'binomial');
        
        % calculate predicted value for logistic model with coefficients b_coeff at
        % the points defined in predVar (for the predicted model)
        yfit = glmval(b_coeff,predVar,'logit');
        
        % yfit is weight associated with model in second column
        % first column weight is 1-yfit
        modelWeights = [1-yfit yfit];
        
    case 'continuous'
        
        % calculate model weights according to AIC values
        minVal = min(AICs,[],2);
        dAIC = AICs - repmat(minVal, 1, size(AICs,2));
        expdAIC = exp(-0.5.*dAIC);
        AIC_weights = expdAIC ./ repmat(sum(expdAIC,2),1,size(expdAIC,2));
        
        % fit logistic model to these weights
        b_coeff = glmfit(trainVar, [AIC_weights(:,2) ones(size(AIC_weights(:,2)))],'binomial');
        
        % calculate predicted value for logistic model with coefficients b_coeff at
        % the points defined in predVar (for the predicted model)
        pred_weights = glmval(b_coeff,predVar,'logit');
        
        modelWeights(:,2) = pred_weights;
        modelWeights(:,1) = 1 - modelWeights(:,2);
        
        
end
