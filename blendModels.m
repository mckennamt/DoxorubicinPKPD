

function [blended, blended_smooth] = blendModels(modelWeights, ...
    modelPredictions, modelPredictions_smooth)


blended = cell(size(modelPredictions,1),1);
blended_smooth = cell(size(modelPredictions,1),1);

if size(modelWeights,2) == 2
    
    % go through each condition in the test set
    for testingIter = 1:size(blended,1)
        
        blended{testingIter} = ...
            modelPredictions{testingIter,1}.*modelWeights(testingIter,1) + ...
            modelPredictions{testingIter,2}.*modelWeights(testingIter,2);
        
        blended_smooth{testingIter} = ...
            modelPredictions_smooth{testingIter,1}.*modelWeights(testingIter,1) + ...
            modelPredictions_smooth{testingIter,2}.*modelWeights(testingIter,2);
        
    end
    
end