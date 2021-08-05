
function Y = runGrowthModel_allData(parms,cellData)

% Optimize this function to get estimates of kp and theta (carrying
% capacity) using pre-treatment and control data
% Uses data from all exposure time experiments

%% Forward solve for pre-treatment data in all wells

Y_est = [];
Y_meas = [];
Y_weight = [];

for eti = 1:length(cellData.ExposureTimes)
    
    %cellData.allYd_tv
    treat_time = cellData.DrugAdded_Tp;
    preTreatData = cellfun(@(x) x(1:treat_time,:),...
        cellData.allYall{eti},'UniformOutput',false);
    preTreatData = cat(2,preTreatData{:});
    Y_meas = [Y_meas; preTreatData(:)];
    
    Tspan = cellData.Tspan{eti};
    Tspan = Tspan(1:treat_time);
    
    Y_cell = cell(size(preTreatData,2),1);
    
    S0_treat = preTreatData(1,:);
    
    for wellIter = 1:size(preTreatData,2)
        %for wellIter = 1:length(treatmentWells)
        
        %[~, Ypr] = ode23(@(t,y) growthModel(t,y,parms),...
        %    Tspan,S0_treat(wellIter));
        Ypr = growthModel_analytic(Tspan, S0_treat(wellIter), parms);
        Y_cell{wellIter} = Ypr;
        
    end
    
    Y_est1 = zeros(size(preTreatData,2)*length(Tspan),1);
    for fillIter = 1:size(preTreatData,2)
        Y_est1(length(Tspan)*(fillIter-1)+1:length(Tspan)*fillIter) = Y_cell{fillIter};
    end
    Y_est = [Y_est;Y_est1];
    Y_weight = [Y_weight; ones(size(Y_est1))];
    
    %% Forward solve for entire timecourse in control wells (no drug)
    % only looking from first post-treatment timepoint onward (drop post-drug
    % even in control wells artificially decreases growth rate)
    
    startTP = cellData.DrugAdded_Tp + 1;
    endTP = cellData.goodControlData;
    
    controlData = cellfun(@(x) x(startTP:endTP,:),...
        cellData.allYall{eti}(1),'UniformOutput',false);
    controlData = controlData{1};
    Y_meas = [Y_meas; controlData(:)];
    
    Tspan = cellData.Tspan{eti};
    Tspan2 = Tspan(startTP:endTP);
    
    Y_cell = cell(size(controlData,2),1);
    Y_weights_cell = cell(size(controlData,2),1);
    
    S0_control = controlData(1,:);
    
    for wellIter = 1:size(controlData,2)
        %for wellIter = 1:length(controlWells)
        
        %[~,Ypr] = ode23(@(t,y) growthModel(t,y,parms),...
        %    Tspan2,S0_control(wellIter));
        Ypr = growthModel_analytic(Tspan2, S0_control(wellIter), parms);
        Y_cell{wellIter} = Ypr;
        
        c_weights = ones(size(Ypr));
        if isfield(cellData, 'allYweights')
            lvtp = cellData.allYweights{eti}{1}(wellIter);
            if lvtp < length(Tspan)
                c_weights = zeros(size(Ypr));
                c_weights(1:(lvtp-cellData.DrugAdded_Tp)) = 1;
            end
        end
        Y_weights_cell{wellIter} = c_weights;
        
    end
    
    Y_est2 = zeros(size(controlData,2)*length(Tspan2),1);
    Y_weight2 = zeros(size(controlData,2)*length(Tspan2),1);
    for fillIter = 1:size(controlData,2)
        Y_est2(length(Tspan2)*(fillIter-1)+1:length(Tspan2)*fillIter) = Y_cell{fillIter};
        Y_weight2(length(Tspan2)*(fillIter-1)+1:length(Tspan2)*fillIter) = Y_weights_cell{fillIter};
    end
    Y_est = [Y_est;Y_est2];
    Y_weight = [Y_weight; Y_weight2];
    
end

Y = (Y_est - Y_meas)./Y_meas .* Y_weight;

end

%%
function dSdt = growthModel(t,S,parms)

% extract parameters
% cell population parameters
kp = parms(1);
theta = parms(2);

% cell populations
dSdt = kp*S(1)*(1-S(1)/theta);

dSdt = dSdt';
end

function N = growthModel_analytic(t,N0,parms)

% extract parameters
% cell population parameters
kp = parms(1);
theta = parms(2);

t = (t-min(t));

N = theta*N0 ./ (N0 + (theta-N0).*exp(-kp.*t));

end