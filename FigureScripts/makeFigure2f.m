close all;clear all;clc


allCellLines = {'SUM149','MDAMB231','MDAMB453','MDAMB468'};

estim =     [.01 .02 10];
lowBound = [0 0 0];
highBound = [Inf Inf Inf];

hold_fits = zeros(3,3,length(allCellLines));
R2 = zeros(1,length(allCellLines));
rMSE = zeros(1,length(allCellLines));

for cellLineIter = 4%1:length(allCellLines)
    
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
    parmEsts{3} = linspace(.1, 1, npe);
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
    
    
    %mean square error
    v= (length(Res)-length(parm));
    rMSE(cellLineIter) = sqrt(sum(Rn)/v);
    
    
    % calculate confidence intervals for each parameter
    dof= (length(Res)-length(parm));
    v = Rn/dof;
    % the parameter covariance matrix
    C = v*inv(J'*J);
    
    % extract the standard errors for each fitted parameter
    se_g = sqrt(diag(C));
    %if cellLineIter==2
    %    se_g = sqrt(abs(diag(C)));
    %else
    %end
    
    % compute the 95% confidence intervals
    tst = tinv(0.975,dof);
    %disp(['95% confidence intervals of estimates of parameters'])
    
    % Normalize by Volume
    kEF = parm(1) * v_I/v_E;
    kFE = parm(1)/parm(2);
    kFB = parm(3);
    
    z_ci = zeros(3,3);
    z_ci(:,1) = [kEF - tst*se_g(1) * v_I/v_E;
        (parm(1)- tst*se_g(1))/(parm(2) + tst*se_g(2));%*parm(1)/parm(2);
        kFB - tst*se_g(3)];
    z_ci(:,2) = [kEF;kFE;kFB];
    z_ci(:,3) = [kEF + tst*se_g(1) * v_I/v_E;
        (parm(1)+ tst*se_g(1))/(parm(2) - tst*se_g(2));%*parm(1)/parm(2);
        kFB + tst*se_g(3)];
    %z_ci = [parm'-tst*se_g,parm',parm'+tst*se_g];
    
    hold_fits(:,:,cellLineIter) = z_ci;
    
end
%% plot it

cf = figure(2);
clf

fs = 12;
min_y = 0;
sf = 1e-5;
yal = [0 2.5*sf;
    0 1.8;
    0 0.3];

allCellLines = {'SUM-149PT','MDA-MB-231','MDA-MB-453','MDA-MB-468'};
cellMarks = {'rx','bo','g^','ks'};

for cellLineIter = 1:length(allCellLines)
    % plot values
    for spi = 1:3
        subplot(1,3,spi)
        %semilogy(-10, 1e10,'k.')
        
        axis([.4 4.6 yal(spi,1) yal(spi,2)])
        
        mv = hold_fits(spi,2, cellLineIter);
        lb = hold_fits(spi,2, cellLineIter) - hold_fits(spi,1, cellLineIter);
        ub = hold_fits(spi,3, cellLineIter) - hold_fits(spi,2, cellLineIter);
        %lb(lb<=min_y) = min_y;
        
        hold on
        ofst = .5;
        cAx = gca;
        h = patch([cellLineIter-ofst cellLineIter+ofst, ...
            cellLineIter+ofst cellLineIter-ofst],...
            [min_y min_y mv mv],[0 0 0],'FaceAlpha',.1,'LineWidth',1,...
            'EdgeColor',[0 0 0]);
        
        errorbar(cellLineIter, mv, lb, ub, 'k',...
            'MarkerSize',10, 'LineWidth',2)
    end
end

for spi = 1:3
    subplot(1,3,spi)
    ax = gca;
    set(gca,'XTick',1:length(allCellLines),'XTickLabel',allCellLines,...
        'XTickLabelRotation',15,'FontSize',fs)
    
    if spi == 1
        pause(.0001)
        set(ax,'YTickLabel',sprintf('%2.1f\n',ax.YTick./sf))
        text(ax.XLim(1), ax.YLim(2)*1.04, 'x10^{-5}','FontSize',fs);
    
    elseif spi==2
        pause(.0001)
        set(ax,'YTickLabel',sprintf('%2.1f\n',ax.YTick))
    
    elseif spi == 3
        pause(.0001)
        set(ax,'YTickLabel',sprintf('%2.1f\n',ax.YTick./.1))
        text(ax.XLim(1), ax.YLim(2)*1.04, 'x10^{-1}','FontSize',fs);
    end
end

set(cf,'Position',[1430 300 1500 350])
%%
R2