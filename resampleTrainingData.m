function [resampd_TrainingData,resampd_WeightData] =...
    resampleTrainingData(allTrainingData,allWeights)

% resample training data for bootstrap analysis

resampd_TrainingData = cell(size(allTrainingData));
resampd_WeightData = cell(size(allWeights));

for treatConds = 1:size(allTrainingData,1)
    
    cData = allTrainingData{treatConds};
    
    % six samples from each, along each column
    [resampd_TrainingData{treatConds},sampdInds] = ...
        datasample(cData,6,2,'Replace',true);
    
    if ~isempty(allWeights)
        cWeights = allWeights{treatConds};
        resampd_WeightData{treatConds} = cWeights(sampdInds);
    end
    
end

if isempty(allWeights)
    resampd_WeightData = [];
end