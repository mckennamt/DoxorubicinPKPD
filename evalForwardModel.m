

function [allPredictions_cell, allPredictions_smooth_cell,Tspan_mod] = ...
    evalForwardModel(predictedParms,...
    fixedParms, allTestData, allTspan_test, postTreat_TP)

allPredictions_cell = cell(size(allTestData));
allPredictions_smooth_cell = cell(size(allTestData));
Tspan_smooth_cell = cell(size(allTestData));

parfor treatCondIter = 1:size(predictedParms,1)
    
    cData = allTestData{treatCondIter};
    Tspan = allTspan_test{treatCondIter};
    Tspan = Tspan(postTreat_TP:end);
    Tspan_mod = linspace(Tspan(1),Tspan(end),100);
    
    fitParms = predictedParms(treatCondIter,:);
    
    Yd = [];
    time_vector = [];
    maxDrugConc = [];
    
    allPredictions = zeros(length(Tspan), size(cData,2));
    allPredictions_smooth = zeros(length(Tspan_mod), size(cData,2));
    
    for initConditIter = 1:size(cData,2)
        
        S0 = cData(postTreat_TP,initConditIter);
        
        % same timepoints as measured data
        %[Tpr, Ypr] = ode45(@(t,y) damageModel_Simplified_TimeMod(t,y,...
        %    fitParms,fixedParms,maxDrugConc,Yd,time_vector),Tspan,S0);
        Ypr = damageModel_analytic(S0, Tspan, fitParms, fixedParms);
        allPredictions(:,initConditIter) = Ypr;
        
        
        %[Tpr_smooth, Ypr_smooth] = ode45(@(t,y) damageModel_Simplified_TimeMod(t,y,...
        %    fitParms,fixedParms,maxDrugConc,Yd,time_vector),Tspan_mod,S0);
        Ypr_smooth = damageModel_analytic(S0, Tspan_mod, fitParms, fixedParms);
        allPredictions_smooth(:,initConditIter) = Ypr_smooth;
         
    end
    
    allPredictions_cell{treatCondIter} = allPredictions;
    allPredictions_smooth_cell{treatCondIter} = allPredictions_smooth;
    Tspan_smooth_cell{treatCondIter} = Tspan_mod;
    
end

Tspan = allTspan_test{1};
Tspan = Tspan(postTreat_TP:end);
Tspan_mod = linspace(Tspan(1),Tspan(end),100);