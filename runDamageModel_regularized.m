
function [Y,J] = runDamageModel_regularized(parms, currentData, allTspan_training, ...
    fixedParms, postTreat_TP, regConstant, regularizer, calcJacobian, currentWeights)

% currentData is cell vector with length corresponding to all
% concentrations at a single exposure time, each cell contains all
% replicates at that given concentration*exposure time

Yd = [];
time_vector = [];
maxDrugConc = [];

% resize parameters into number of concentrations * number parameters
% matrix
np = length(parms)/length(currentData);
inParms = reshape(parms,length(currentData),np);

dataVector = cell(size(currentData));
modelEstVector = cell(size(currentData));
weightVector = cell(size(currentData));
JVector = cell(size(currentData));

numTreatConds = length(currentData);

% run each concentration in dataset forward
for treatmentIter = 1:length(currentData)
    
    cYall = currentData{treatmentIter};
    cYall = cYall(postTreat_TP:end,:);
    
    Tspan = allTspan_training{treatmentIter};
    Tspan_modified = Tspan(postTreat_TP:end);
    
    if isempty(currentWeights)
        cWeights = length(Tspan_modified) .* ones(1, size(cYall,2));
    else
        cWeights = currentWeights{treatmentIter};
        cWeights = cWeights - (postTreat_TP-1);
    end
    
    c_parms = inParms(treatmentIter,:);
    
    [dataVector{treatmentIter}, modelEstVector{treatmentIter}, ...
        weightVector{treatmentIter}, Jcalc] = ...
        forwardSolve_DamageModel_Simplified_TimeMod(c_parms,...
        cYall, Tspan_modified, fixedParms, maxDrugConc, Yd, time_vector,calcJacobian,cWeights);
    
    if calcJacobian
        numSamples_perConc = numel(cYall);
        sr = (treatmentIter-1)*numSamples_perConc + 1;
        er =  treatmentIter*numSamples_perConc;
        allR = repmat([sr:er]',np,1);
        allC = repmat(treatmentIter:numTreatConds:length(parms),length(sr:er),1);
        JVector{treatmentIter} = [allR allC(:) Jcalc(:)];
    end
    
    
    
end

% construct vector with all data, model predictions, and weights
Ydata = cat(1,dataVector{:});
Yest = cat(1,modelEstVector{:});
Yweight = cat(1,weightVector{:});
%Yweight(:) = 1;

% calculate smoothness of parameter maps

% tikhonov regularization ||Ydata-Yest||2 + alpha*||T*parms||
% enforce smoothness on data
regFunVal = regulizerFun(inParms,regularizer);

Y = [(Ydata - Yest).*Yweight./Ydata; ...
    regConstant.*regFunVal];

if calcJacobian
    cYall = currentData{1};
    cYall = cYall(postTreat_TP:end,:);
    numSamples_perConc = numel(cYall);
    regRows =  repmat(length(currentData)*numSamples_perConc + [1:length(parms)]', length(parms), 1);
    regCols =  repmat([1:length(parms)]', length(parms), 1);
    regJ = repmat(regConstant./regularizer(:,2), length(parms), 1);
    
    J = cat(1,JVector{:});
    J = [J; regRows regCols regJ];
    J = sparse(J(:,1), J(:,2), J(:,3), max(J(:,1)), max(J(:,2)));
else
    J = [];
end

end


function [Ydata, Yest, Yweight,Jacob] = forwardSolve_DamageModel_Simplified_TimeMod(parms,Yall,Tspan,...
    fixedParms,maxDrugConc,Yd,time_vector, calcJacobian, cWeights)

% Solve ODE
%Tspan = cellData.times(cellData.DrugAdded_Tp:end);
%Yall is matrix length(Tspan)xNumberReplicates

Ymod = cell(1,size(Yall,2));
Jmod = cell(1,size(Yall,2));

for wellIter = 1:size(Yall,2)
    
    S0 = Yall(1,wellIter);
    
    %[~, Ycell] = ode45(@(t,y) damageModel_Simplified_TimeMod(t,y,...
    %    parms,fixedParms,maxDrugConc,Yd,time_vector),...
    %    Tspan,S0);%ode23
    if calcJacobian
        [Ycell,Jcell] = damageModel_analytic(S0, Tspan, parms, fixedParms);
        Jmod{wellIter} = Jcell;
    else
        Ycell = damageModel_analytic(S0, Tspan, parms, fixedParms);
    end
    
    Ymod{wellIter} = Ycell;
end

Jacob = cat(1, Jmod{:});

Yest = zeros(size(Yall));
Yweight = zeros(size(Yall));
goodData = true(size(Yall));
lTspan = length(Tspan);
%Y = zeros(6*length(Tspan),1);
for wellIter = 1:size(Yall,2)
    Yest(:,wellIter) = Ymod{wellIter};
    
    %weightVec = ones(size(Yall,1),1);
    weightVec = zeros(size(Yall,1),1);
    weightVec(1:cWeights(wellIter)) = 1;
    
    %reachCC = Yall(:,wellIter) >= .8*fixedParms(2);
    %fcc = find(reachCC,1);
    %lcc = find(reachCC,1,'last');
    %if ~isempty(lcc) && ~isempty(fcc) && lcc<lTspan
    %    if lcc < fcc+3 && fcc+3<lTspan
    %        %goodData(fcc+4:end,wellIter) = false;
    %        goodData(fcc+4:end,wellIter) = true;
    %    else
    %        %goodData(lcc+1:end,wellIter) = false;
    %        goodData(fcc+4:end,wellIter) = true;
    %    end
    %end
    
    %weightVec(8) = 0;
    Yweight(:,wellIter) = weightVec;
    
end

Yest(Yest<=0) = 1;

Yest = Yest(:);
Ydata = Yall(:);
Yweight = Yweight(:);

Yest = Yest(goodData(:));
Ydata = Ydata(goodData(:));
Yweight = Yweight(goodData(:));

%Y = (Ydata - Yest).*Yweight./Ydata;

end

function regFunVal = regulizerFun(inParms,regularizer)

regFunVal = zeros(2*numel(inParms),1);

switch regularizer{1}
    
    case 'small'
        hld_normalizer = regularizer{2};
        
        regFunVal(1:numel(inParms)) =...
            (inParms(:)-hld_normalizer(:,1))./hld_normalizer(:,2);
        
    case 'smooth'
        % ensure it's smooth
        hld_normalizer = regularizer{2};
        doxConc = hld_normalizer(1:size(inParms,1),1);
        
        [doxConc, srtInds] = sort(doxConc);
        
        normalizeParms = regularizer{3};
        inParms = inParms./repmat(normalizeParms, length(doxConc), 1);
        
        inParms = inParms(srtInds,:);
        
        
        fEnt = inParms(1,:) - (inParms(2,:)-inParms(1,:));
        lEnt = inParms(end,:) + (inParms(end,:)-inParms(end-1,:));
        inParms_m1 = [fEnt; inParms(1:end-1,:)];
        inParms_p1 = [inParms(2:end,:); lEnt];
        
        h = diff(doxConc);
        h = [h(1); h] + [h; h(end)];
        h = repmat(h, 1, size(inParms,2));
        
        hld_regFunVal(srtInds,:) = (inParms_p1-inParms_m1)./(2.*h);
        
        %if size(inParms,2) == 2
        %    hld_regFunVal(srtInds,:) = (inParms_p1-inParms_m1)./(2.*h);
        %else
        %    hld_regFunVal(srtInds,:) = (inParms_p1-inParms_m1)./(2.*h) + ...
        %        [abs(inParms(:,1)), .01./abs(inParms(:,2)), zeros(size(inParms,1), 1)];
        %end
        
        
        regFunVal(1:numel(inParms)) = hld_regFunVal(:);
        
        if size(hld_normalizer,2) == 2
            doxConc = hld_normalizer(1:size(inParms,1),2);
            
            [doxConc, srtInds] = sort(doxConc);
            
            h = diff(doxConc);
            h = [h(1); h] + [h; h(end)];
            h = repmat(h, 1, size(inParms,2));
            
            hld_regFunVal(srtInds,:) = (inParms_p1-inParms_m1)./(2.*h);
            %if size(inParms,2) == 2
            %    hld_regFunVal(srtInds,:) = (inParms_p1-inParms_m1)./(2.*h);
            %else
            %    hld_regFunVal(srtInds,:) = (inParms_p1-inParms_m1)./(2.*h) + ...
            %        [abs(inParms(:,1)), .01./abs(inParms(:,2)), zeros(size(inParms,1), 1)];
            %end
            
            regFunVal((numel(inParms)+1):end) = hld_regFunVal(:);
            
        end
        
end

end