close all;%clear all;clc
clearvars -except AIC_compartmentModel


allCellLines = {'SUM149','MDAMB231','MDAMB453','MDAMB468'};

estim =     [.01 .02];% 10 10];
lowBound = [0 0];% 0 0];
highBound = [Inf Inf];% Inf Inf];

hold_fits = zeros(3,3,length(allCellLines));
R2 = zeros(1,length(allCellLines));
rMSE = zeros(1,length(allCellLines));



for cellLineIter = 1:length(allCellLines)
    
    % load current cell line data
    CellLine = allCellLines{cellLineIter};
    %load([CellLine 'CompartmentParameters.mat'])
    if cellLineIter==1
        load(['IntensityTimecourses_20161022_' CellLine '.mat']);
    else
        load(['IntensityTimecourses_' CellLine '.mat']);
    end
    
    v_E = 250; %uL
    % cell volume (um^3) * number cells * 1000^-3 (um3 to mm3)
    % 1mm3 = 1uL
    v_I = getCellVolume(CellLine) * 10000 * 1000^-3;%uL
    
    % clean up data
    allTimes_hrs = allTimes_hrs(1:end-1);
    concMatrix = concMatrix(:,1:end-1,:);
    allTimes_hrs_mod = allTimes_hrs + abs(min(allTimes_hrs));
    analysisWells = [1:8 11:18];
    
    % get range of initial estimates and set bounds on parameter estimates
    npe = 3;
    clear parmEsts
    parmEsts{1} = linspace(.01, .2, npe);
    parmEsts{2} = linspace(.01, .2, npe);
    %parmEsts{3} = linspace(.1, 1, npe);
    %parmEsts{4} = linspace(.01, 1, npe);
    hldGrid = cell(1,numel(parmEsts));
    [hldGrid{:}] = ndgrid(parmEsts{:});
    hldGrid = cellfun(@(x) x(:), hldGrid,'UniformOutput',false);
    startPointList = cat(2,hldGrid{:});
    custpts = CustomStartPointSet(startPointList);
    
    % set up optimization parameters
    lsqOpts = optimoptions('lsqnonlin','MaxFunEvals',500,'Display','none',...
        'TolFun', 1e-6,'MaxIter', 5000,'TolX', 1e-5);
    
    % set up multistart optimization problem
    prob = createOptimProblem('lsqnonlin','x0',startPointList(1,:),...
        'objective',@(kpars) runCompartmentModel_CellLine(kpars,concMatrix,...
        allTimes_hrs_mod, analysisWells, aifDuration,startDrug, endDrug),...
        'lb',lowBound,'ub',highBound,'options',lsqOpts);
    ms = MultiStart('UseParallel',true,'Display','iter');
    
    % run optimization
    [parm,~, ~, ~,ms_solutions] = run(ms,prob,custpts);
    
    % run at the solution to calculate the residual and jacobian to calculate
    % the confidence intervals for each parameter
    [parm,Rn,Res,whyExit,~,~,J] = lsqnonlin(@(kpars) runCompartmentModel_CellLine(kpars,concMatrix,...
        allTimes_hrs_mod, analysisWells, aifDuration,startDrug, endDrug),...
        parm, lowBound, highBound, lsqOpts);
    
    %R2(cellLineIter) = runCompartmentModel_r2(parm,concMatrix,...
    %    allTimes_hrs_mod, analysisWells, aifDuration,startDrug,endDrug);
    [~,R2(cellLineIter)] = runCompartmentModel_CellLine(parm,concMatrix,...
        allTimes_hrs_mod, analysisWells, aifDuration,startDrug,endDrug);
    
   
    [Y_meas,Y_pred,Y_weight] = runCompartmentModel_CellLine_forResidCalc(parm,concMatrix,...
        allTimes_hrs_mod, analysisWells, aifDuration,startDrug,endDrug);
    
    n = sum(sum(Y_meas>0));
    resNorm = sum((Y_meas(~isnan(Y_meas)) - Y_pred(~isnan(Y_meas))).^2.*Y_weight(~isnan(Y_meas)));
    k = 3;
    AIC_compartmentModel(cellLineIter,3) = n*log(resNorm/n)+(2*k)+(((2*k)*(k+1))/(n-k-1));
    
end
