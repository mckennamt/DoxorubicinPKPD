clear all;clc

load('DataToRuneFAST_prevFig5_OptRegConst_20161213.mat')

rng('default');
addpath('./Data')

poolobj = gcp('nocreate'); % If no pool, do not create new one.
if isempty(poolobj)
    c = parcluster('local');
    poolsize = c.NumWorkers;
    myPool = parpool(poolsize);
else
    myPool = poolobj;
end

addpath('./SensitivityAnalysis')
addpath('./SensitivityAnalysis/eFAST')
addpath('./SensitivityAnalysis/GSAT')
addpath('./SensitivityAnalysis/lhs-prcc')

%%
close all

disp('Ready to go!')

pause(1)

modNames = {'5A','5B'};
for cellLineIter = 2%[2 3 1 4]
    for modIter = 1:2
        
        cellData = allCellData{cellLineIter};
        cfits = collect_Fits{cellLineIter};
        c_AIC = collect_AIC{cellLineIter};
        
        %[holdSens_r, holdSens_kdb, holdSens_kda, holdSens_wb,l_cb_smoove] = ...
        %    eFAST_DoseResponse(cellData,cfits,c_AIC,allCellLine{cellLineIter});
        
        modelSwitch = modNames{modIter};
        eFAST_DoseResponse_Take2(cellData,cfits,allCellLine{cellLineIter},modelSwitch);
        
    end
end

%eFAST_DoseResponse_SupplementMaterials('General');