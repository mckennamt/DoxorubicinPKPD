

% use other cell lines to set fitting parameters

function fitSettings = optimizeFittingParameters_only24(otherCellLines, funToFit)

% get parameter values for all cell lines
numberCellLines = length(otherCellLines);

% regularization constant for each model, and span for all models (only
% want to identify a single set of training data)
numFittingParms = length(funToFit) + 1;

% each parameter can use anywhere from 10%-100% of the data in the local
% regression
allSpan = [.2:.05:.9]';
allReg = 1:10;%:2:20

[fitParmMat{1:numFittingParms}] = ndgrid(allReg, allReg, allSpan);

% get matrix where each row is unique set of span combinations
allFittingParmCombos = (cellfun(@(x) x(:), fitParmMat,'UniformOutput',false));
allFittingParmCombos = cat(2,allFittingParmCombos{:});

%% loop through all span combinations to find best fitting parameters

allSpanErrors = zeros(size(allFittingParmCombos,1), length(otherCellLines));

uregCombos = unique(allFittingParmCombos(:,1:length(funToFit)),'rows');

for cellLineIter = 1:numberCellLines
    
    CellLine = otherCellLines{cellLineIter};
    fprintf('%s\n', CellLine)
    for comboIter = 1:size(uregCombos,1)
        fprintf('%d ', comboIter)
        c_fittingParms = uregCombos(comboIter,:);
        
        [cellData, cell_AIC,...
            cellFits, ~, cell_kpTheta, ~] = ...
            fitChemoModel_cellLine(CellLine, funToFit, c_fittingParms);
        %keyboard
        
        spans = allFittingParmCombos(ismember(allFittingParmCombos(:,1:length(funToFit)),...
            c_fittingParms,'rows'),3);
        
        hold_allSpanErrors = zeros(size(spans));
        
        parfor span_iter = 1:length(spans)
        %for span_iter = 1:length(spans)
            c_span = spans(span_iter);
            
            
            MPE_cl = zeros(1,2);
            ioi_tests = logical([1 0 0;0 1 0]');
            
            for testIter = 1:2
                
                % collect Training Data
                ioi_train = logical([0 0 1]');
                
                alltrainParms = cellFits(ioi_train,:);
                alltrainAICs = cell_AIC(ioi_train,:);
                trainParms = cell(1, length(funToFit));
                trainAICs = zeros(sum(ioi_train)*size(cell_AIC{1},1), size(cell_AIC,2));
                for catIter = 1:length(funToFit)
                    trainParms{catIter} = cat(1, alltrainParms{:,catIter});
                    trainAICs(:,catIter) = cat(1, alltrainAICs{:,catIter});
                end
                
                % remove control data (no drug added)
                cb_training = cat(2,cellData.dox_pss{ioi_train});
                cb_training = cb_training(2:end,:);
                cb_training = cb_training(:);
                
                auc_training = cat(2,cellData.auc{ioi_train});
                auc_training = auc_training(2:end,:);
                auc_training = auc_training(:);
                
                predictorVariables_training = [cb_training auc_training];
                
                % sort training data (for regularization calculation)
                [predictorVariables_training, srtedInds] = sortrows(predictorVariables_training);
                trainParms = cellfun(@(x) x(srtedInds,:), trainParms,'UniformOutput',false);
                trainAICs = trainAICs(srtedInds,:);
                
                % collect Test Data
                ioi_test = ioi_tests(:,testIter);
                
                postTreat_TP = cellData.DrugAdded_Tp + 1;
                
                % remove control data (no drug added)
                allTestData = cat(2,cellData.allYall{ioi_test});
                allTestData = allTestData(2:end,:);
                allTestData = allTestData(:);
                %allTestData = cellfun(@(x) x(postTreat_TP:end,:), allTestData,...
                %    'UniformOutput',false);
                
                if isfield(cellData, 'allYweights')
                    allTestWeights = cat(2,cellData.allYweights{ioi_test});
                    allTestWeights = allTestWeights(2:end,:);
                    allTestWeights = allTestWeights(:);
                else
                    allTestWeights = [];
                end
                
                cb_test = cat(2,cellData.dox_pss{ioi_test});
                cb_test = cb_test(2:end,:);
                cb_test = cb_test(:);
                
                auc_test = cat(2,cellData.auc{ioi_test});
                auc_test = auc_test(2:end,:);
                auc_test = auc_test(:);
                
                predictorVariables_test = [cb_test auc_test];
                
                % define timespan for model
                allTspan_test = cellData.Tspan(ioi_test);
                allTspan_test = repmat(allTspan_test,size(allTestData,1),1);
                allTspan_test = allTspan_test(:);
                %allTspan_test = cellfun(@(x) x(postTreat_TP:end), allTspan_test,...
                %    'UniformOutput',false);
                
                if isempty(allTestWeights)
                    weightCell = cellfun(@(x) ones(size(x)), allTestData, 'UniformOutput',false);
                else
                    weightCell = cellfun(@(x) zeros(size(x)), allTestData, 'UniformOutput',false);
                    for cellIter = 1:length(weightCell)
                        for colIter = 1:size(weightCell{cellIter},2)
                            weightCell{cellIter}(1:allTestWeights{cellIter}(colIter),colIter) = 1;
                        end
                    end
                end
                
                % calculate mean percent error for each test set
                MPE_cl(testIter) = calculate_MeanPercentError(...
                    trainParms, trainAICs, predictorVariables_training,...
                    cell_kpTheta, funToFit, c_span,...
                    allTestData, allTspan_test, predictorVariables_test,...
                    postTreat_TP,weightCell);
                
            end
            
            %ioi = ismember(allFittingParmCombos, [c_fittingParms c_span],'rows');
            %allSpanErrors(ioi,cellLineIter) = mean(MPE_cl);
            hold_allSpanErrors(span_iter) = mean(MPE_cl);
        end
        ioi = ismember(allFittingParmCombos, [repmat(c_fittingParms, length(spans), 1) spans],'rows');
        allSpanErrors(ioi,cellLineIter) = hold_allSpanErrors;
    end
    fprintf('\n')
end

[~,minInd] = min(mean(allSpanErrors,2));
fitSettings = allFittingParmCombos(minInd,:);