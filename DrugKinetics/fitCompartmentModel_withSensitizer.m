
clear all;clc
CellLine = 'SKMEL5';
currDate = 20161110;

%load(['IntensityTimecourses_' CellLine '.mat']);
%load(['IntensityTimecourses_' CellLine '_RFPTest.mat']);
%load(['IntensityTimecourses_20161022_' CellLine '.mat']);
load(['IntensityTimecourses_' num2str(currDate) '_' CellLine '.mat'])

v_E = 250; %uL
% cell volume (um^3) * number cells * 1000^-3 (um3 to mm3)
% 1mm3 = 1uL
v_I = getCellVolume('MDAMB468') * 10000 * 1000^-3;%uL


% fit compartment model

allTimes_hrs = allTimes_hrs(1:end-1);
concMatrix = concMatrix(:,1:end-1,:);


allTimes_hrs_mod = allTimes_hrs + abs(min(allTimes_hrs));
%options=optimset('TolFun',1e-6,'TolX',1e-6,'MaxIter',1e50,'Display','off','PrecondBandWidth',0);

analysisWells = [11:19];
%analysisWells = [11:18]% 11:18];

lowBound = [0 0 0];
%highBound = [10 10 10 10 1e5];
highBound = [Inf Inf Inf];

% 3 compartment 3 parameter
% kpars =   [k1 k2 k3]
%estim =     [.01 .02 10 .05];
%estim =     [2 0.05 .05];
%lowBound =  [0.001  0.001 .001];
%highBound = [10  5 5];

npe = 5;
% create parameter estimates
clear parmEsts
parmEsts{1} = linspace(.1, .5, npe);%5
parmEsts{2} = linspace(.1, .5, npe);%5
parmEsts{3} = linspace(.01, .2, npe);%5
%parmEsts{4} = linspace(50, 100, npe);%5
%parmEsts{5} = linspace(1, 10, npe);%5
hldGrid = cell(1,numel(parmEsts));
[hldGrid{:}] = ndgrid(parmEsts{:});
hldGrid = cellfun(@(x) x(:), hldGrid,'UniformOutput',false);
startPointList = cat(2,hldGrid{:});
custpts = CustomStartPointSet(startPointList);

% set up multi-start optimization
lsqOpts = optimoptions('lsqnonlin','MaxFunEvals',100,'Display','none',...
    'TolFun', 1e-6,'MaxIter', 5000,'TolX', 1e-8);
prob = createOptimProblem('lsqnonlin','x0',startPointList(1,:),...
    'objective',@(kpars) runCompartmentModel_CellLine(kpars,concMatrix,...
    allTimes_hrs_mod, analysisWells, aifDuration,startDrug, endDrug),...
    'lb',lowBound,'ub',highBound,'options',lsqOpts);
ms = MultiStart('UseParallel',true,'Display','iter');
[parm,~, ~, ~,ms_solutions] = run(ms,prob,custpts);
disp(parm)
% run at the solution to calculate the residual and jacobian to calculate
% the confidence intervals for each parameter
[parm,Rn,Res,whyExit,~,~,J] = lsqnonlin(@(kpars) runCompartmentModel_CellLine(kpars,concMatrix,...
    allTimes_hrs_mod, analysisWells, aifDuration,startDrug, endDrug),...
    parm, lowBound, highBound, lsqOpts);

[~,R2] = runCompartmentModel_CellLine(parm,concMatrix,...
        allTimes_hrs_mod, analysisWells, aifDuration,startDrug,endDrug)

% Calculate the confidence interval for each parameter
dof= (length(Res)-length(parm));
v = Rn/dof;
C = v*inv(J'*J); % the parameter covariance matrix
se_g = sqrt(diag(C)); % extract the standard errors for each fitted parameter
tst = tinv(0.975,dof); % compute the 95% confidence intervals
z_ci = [parm'-tst*se_g,parm',parm'+tst*se_g]


kEF = parm(1) * v_I/v_E;
kFE = parm(1)/parm(2);
kFB = parm(3);

%% View compartment model fit for each data well

figure(997);clf
for wellIter = [1:8]%analysisWells
    %for wellIter = [1:8]%analysisWells
    
    
    %     lsqOpts = optimoptions('lsqnonlin','MaxFunEvals',100,'Display','none',...
    %     'TolFun', 1e-50,'MaxIter', 5000,'TolX', 1e-5);
    % prob = createOptimProblem('lsqnonlin','x0',startPointList(1,:),...
    %     'objective',@(kpars) runCompartmentModel_CellLine(kpars,concMatrix,...
    %     allTimes_hrs_mod, wellIter, aifDuration,startDrug, endDrug),...
    %     'lb',lowBound,'ub',highBound,'options',lsqOpts);
    % ms = MultiStart('UseParallel',true,'Display','none');
    % [parm,~, ~, ~,ms_solutions] = run(ms,prob,custpts);
    % disp(wellIter)
    % disp(parm)
    % allParms = [allParms; parm];
    
    cellUptakeSignal = squeeze(concMatrix(wellIter,1:end,1));
    inputFcn = squeeze(concMatrix(wellIter,1:end,2));
    daif = squeeze(concMatrix(wellIter,1:end,3));
    blankFcn = squeeze(concMatrix(wellIter,1:end,4));
    
    if aifDuration(wellIter) == 6
        c_endDrug = endDrug(1);
    elseif aifDuration(wellIter) == 12
        c_endDrug = endDrug(2);
    else
        c_endDrug = endDrug(1);
    end
    
    % clean up uptake curves
    [cellUptakeSignal, inputFcn] = correctUptakeSignal(cellUptakeSignal,...
        inputFcn, blankFcn, startDrug, c_endDrug);
    
    allTimes_upreg = linspace(0,30,100);
    allTimes_hrs_mod = [allTimes_hrs_mod];%; moreTime];
    inputFcn_mod = [inputFcn];% 0.*ones(1,length(moreTime))];
    upregInput = interp1(allTimes_hrs_mod,inputFcn_mod,allTimes_upreg);
    
    % make sure to use same mode as in runCompartmentModel_CellLine
    fwdModEval = modelCt3C3P(parm,upregInput,allTimes_upreg);
    
    %p_mod = parm;
    %p_mod(2) = parm(1)/parm(2);
    %S0 = [0 0];
    %[~, Yfree_Ybound] = ode45(@(t,y) chemoCompartmentModel_odes(t,y,...
    %    p_mod,upregInput,allTimes_upreg),allTimes_upreg,S0);
    %fwdModEval = sum(Yfree_Ybound,2);
    
    %[~,Cout] = ode45(@(t,y) doxCarrierModel(t,y,parm,upregInput, allTimes_upreg),...
    %    allTimes_upreg, [0 0]);
    %fwdModEval = sum(Cout,2);
    
    %[~,fwdModEval] = ode45(@(t,y) twoPeakDoxModel(t,y,parm,upregInput, allTimes_upreg),...
    %    allTimes_upreg, [cellUptakeSignal(2)]);
    
    figure(997);
    if wellIter>10
        wellIter = wellIter - 10;
    end
    subplot(2,4,wellIter)
    hold on
    plot(allTimes_hrs, (cellUptakeSignal), 'ko')
    plot(allTimes_hrs, inputFcn, 'bs')
    plot(allTimes_upreg - allTimes_hrs_mod(startDrug+1), fwdModEval, 'r-')
    plot(allTimes_upreg - allTimes_hrs_mod(startDrug+1), upregInput,'g-')
    legend('Cell uptake','Input Function','Model Fit','Location','Best')
    xlabel('Time (hrs)')
    ylabel('Concentration (nM)')
    hold off
    
    %pause
    
end


%%

%save([CellLine 'CompartmentParameters.mat'],'concMatrix','z_ci',...
%    'allTimes_hrs','allTimes_hrs_mod','analysisWells',...
%    'startDrug','endDrug','aifDuration','drugOnlyConcs',...
%    'allConcs','concFluorInten_Map')

