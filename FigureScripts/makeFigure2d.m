

allCellLines = {'MDAMB468'};

estim =     [.01 .02 10];
lowBound = [0 0 0];
highBound = [1 100 100];

cellLineIter = 1;

% load current cell line data
CellLine = allCellLines{cellLineIter};
load([CellLine 'CompartmentParameters.mat'])

concMatrix(:,end,:) = [];
allTimes_hrs_mod(end) = [];
concFluorInten_Map(end,:) = [];

% get range of initial estimates and set bounds on parameter estimates
parmEsts{1} = linspace(.01, .2, 5);
parmEsts{2} = linspace(.01, .2, 5);
parmEsts{3} = linspace(.1, 20, 5);
hldGrid = cell(1,numel(parmEsts));
[hldGrid{:}] = ndgrid(parmEsts{:});
hldGrid = cellfun(@(x) x(:), hldGrid,'UniformOutput',false);
startPointList = cat(2,hldGrid{:});
custpts = CustomStartPointSet(startPointList);

% set up optimization parameters
lsqOpts = optimoptions('lsqnonlin','MaxFunEvals',10000,'Display','none',...
    'TolFun', 1e-50,'MaxIter', 5000,'TolX', 1e-5);

% limit to just well of interest
concMatrix = concMatrix(1,:,:);
analysisWells = analysisWells(1);
aifDuration = aifDuration(1);
endDrug = endDrug(1);


% set up multistart optimization problem
prob = createOptimProblem('lsqnonlin','x0',startPointList(1,:),...
    'objective',@(kpars) runCompartmentModel_CellLine(kpars,concMatrix,...
    allTimes_hrs_mod, analysisWells, aifDuration,startDrug, endDrug),...
    'lb',lowBound,'ub',highBound,'options',lsqOpts);
ms = MultiStart('UseParallel',true,'Display','off');

% run optimization
[parm,~, ~, ~,ms_solutions] = run(ms,prob,custpts);

% run at the solution to calculate the residual and jacobian to calculate
% the confidence intervals for each parameter
[parm,Rn,Res,whyExit,~,~,J] = lsqnonlin(@(kpars) runCompartmentModel_CellLine(kpars,concMatrix,...
    allTimes_hrs_mod, analysisWells, aifDuration,startDrug, endDrug),...
    parm, lowBound, highBound, lsqOpts);

dof= (length(Res)-length(parm));
v = Rn/dof;
% the parameter covariance matrix
C = v*inv(J'*J);

% extract the standard errors for each fitted parameter

se_g = sqrt(diag(C));

% compute the 95% confidence intervals
tst = tinv(0.975,dof);
%disp(['95% confidence intervals of estimates of parameters'])
z_ci = [parm'-tst*se_g,parm',parm'+tst*se_g];

%%
wellIter = 1;

daif = squeeze(concMatrix(wellIter,1:end,3));
inputFcn = squeeze(concMatrix(wellIter,1:end,2));
blankFcn = squeeze(concMatrix(wellIter,1:end,4));


if aifDuration(wellIter) == 6
    c_endDrug = endDrug(1);
elseif aifDuration(wellIter) == 12
    c_endDrug = endDrug(2);
end

cellUptakeSignal = squeeze(concMatrix(wellIter,1:end,1));

cellUptakeSignal = cellUptakeSignal-blankFcn;
inputFcn = inputFcn-blankFcn;
daif = daif-blankFcn;
cellUptakeSignal(1:end) = cellUptakeSignal(:) - mean(cellUptakeSignal(1:startDrug));

%cellUptakeSignal(1:startDrug) = cellUptakeSignal(1:startDrug) - blankFcn(1:startDrug);
%cellUptakeSignal(c_endDrug+2:end) = cellUptakeSignal(c_endDrug+2:end) - blankFcn(c_endDrug+2:end);
%cellUptakeSignal(c_endDrug+2:end) = cellUptakeSignal(c_endDrug+2:end) - inputFcn(c_endDrug+2:end);

cellUptakeSignal(startDrug+1:c_endDrug+1) = cellUptakeSignal(startDrug+1:c_endDrug+1) - inputFcn(startDrug+1:c_endDrug+1);
%cellUptakeSignal(startDrug+1:c_endDrug+1) = cellUptakeSignal(startDrug+1:c_endDrug+1) - daif(startDrug+1:c_endDrug+1);
%cellUptakeSignal(startDrug+1:c_endDrug+1) = cellUptakeSignal(startDrug+1:c_endDrug+1)+abs(min(cellUptakeSignal));

inputFcn = inputFcn'./concFluorInten_Map(:,1);
inputFcn(1:6) = 0;

allTimes_upreg = linspace(0,allTimes_hrs_mod(end),100);
moreTime = [allTimes_hrs_mod(end)];%+1:1:50]';
allTimes_hrs_mod = [allTimes_hrs_mod];%; moreTime];
inputFcn_mod = [inputFcn];% 0.*ones(1,length(moreTime))];
%inputFcn_mod = zeros(size(inputFcn_mod));
%inputFcn_mod(startDrug:c_endDrug+1) = max(inputFcn(:));
upregInput = interp1(allTimes_hrs_mod,inputFcn_mod,allTimes_upreg);

lb = [z_ci(1,1) z_ci(2,1) z_ci(3,1)];
ub = [z_ci(1,3) z_ci(2,3) z_ci(3,3)];

fwdModEval = modelCt3C3P(parm,upregInput,allTimes_upreg);
fwdModEval_lb = modelCt3C3P(lb,upregInput,allTimes_upreg);
fwdModEval_ub = modelCt3C3P(ub,upregInput,allTimes_upreg);

theShade = zeros(length(allTimes_upreg)*2,2);
theShade(:,1) = [allTimes_upreg'; allTimes_upreg(end:-1:1)'];
theShade(1:length(allTimes_upreg),2) = fwdModEval_lb';
theShade(length(allTimes_upreg)+1:end,2) = fwdModEval_ub(end:-1:1)';

fvad = .8.*ones(size(theShade,1),1);
fvad(length(allTimes_upreg)) = .2;
fvad(end) = .2;

cellUptakeSignal = cellUptakeSignal'./concFluorInten_Map(:,1);
cellUptakeSignal(1:startDrug) = 0;

cf = figure(997);clf
hold on
plot(allTimes_hrs_mod - allTimes_hrs_mod(startDrug), (cellUptakeSignal), 'ko')
plot(allTimes_hrs_mod - allTimes_hrs_mod(startDrug), inputFcn, 'bs-')
%plot(allTimes_upreg - allTimes_hrs_mod(startDrug), fwdModEval, 'r-')

theAlphaVal = .1;
patch(theShade(:,1)-allTimes_hrs_mod(startDrug),theShade(:,2),[1 0 0],'FaceAlpha',theAlphaVal,...
    'EdgeColor',[1 0 0],'EdgeAlpha','flat','LineStyle','-',...
    'FaceVertexAlphaData',fvad)
%plot(allTimes_upreg - allTimes_hrs_mod(startDrug+1), upregInput,'g-')
%xlabel('Time (hrs)')
%ylabel('Concentration (nM)')
hold off
axis([-3 22 -100 2500])
fs = 24;
set(gca,'FontSize',fs)
ax = gca;
set(gca,'YTickLabel',sprintf('%2.1f\n',ax.YTick/1e3))
text(-3, ax.YLim(2)*1.03, 'x10^{3}','FontSize',fs);
%set(cf,'Position',[80 1080 570 425])
set(cf,'Position',[80 1080 800 425])
cl = legend('Intracellular','Extracellular','Model Fit','Location','Best')
set(cl,'FontSize',28)