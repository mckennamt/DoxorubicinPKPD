



function [resNorm,numSamp,percentDiff] = runDamageModel_AICCalc(parms,Yall,Tspan,...
    fixedParms,maxDrugConc,Yd,time_vector,cWeights)

% Solve ODE
%Tspan = cellData.times(cellData.DrugAdded_Tp:end);
%Yall is matrix length(Tspan)xNumberReplicates

Ymod = cell(1,size(Yall,2));
Yweight_cell = cell(1,size(Yall,2));

parfor wellIter = 1:size(Yall,2)
    
    S0 = Yall(1,wellIter);
    %[~, Ycell] = ode45(@(t,y) damageModel_Simplified_TimeMod(t,y,...
    %    parms,fixedParms,maxDrugConc,Yd,time_vector),...
    %    Tspan,S0);%ode23
    
    Ycell = damageModel_analytic(S0, Tspan, parms, fixedParms);
    
    yweight = zeros(size(Ycell));
    yweight(1:cWeights(wellIter)) = 1;
    
    Ymod{wellIter} = Ycell;
    Yweight_cell{wellIter} = yweight;
end

Yest = zeros(size(Yall));
Yweight = zeros(size(Yall));

for wellIter = 1:size(Yall,2)
    Yest(:,wellIter) = Ymod{wellIter};
    Yweight(:,wellIter) = Yweight_cell{wellIter};
end

Yest(Yest<=0) = 1;

Yest = Yest(:);
Ydata = Yall(:);
Yweight = Yweight(:);

%numSamp = length(Yest);
numSamp = sum(Yweight);
resNorm = sum((Ydata - Yest).^2 .* Yweight);
percentDiff_nz = abs(Ydata-Yest)./Ydata;
percentDiff = mean(percentDiff_nz(logical(Yweight)));

