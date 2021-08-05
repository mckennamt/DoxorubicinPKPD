close all;clear all;clc


allCellLines = {'SUM149','MDAMB231','MDAMB453','MDAMB468'};

estim =     [.01 .02 10];
lowBound = [0 0 0];
highBound = [Inf Inf Inf];

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

%%

%close all;
figure(1);clf;
figure(2);clf;
figure(3);clf;
figure(4);clf;

allCellLines = {'SUM149','MDAMB231','MDAMB453','MDAMB468'};

for cellLineIter = 1:length(allCellLines)
    
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
    
    %kEF = parm(1) * v_I/v_E;
    %kFE = parm(1)/parm(2);
    %kFB = parm(3);
    
    parm = [];
    parm(1) = hold_fits(1,2,cellLineIter)/(v_I/v_E);
    parm(2) = parm(1)/hold_fits(2,2,cellLineIter);
    parm(3) = hold_fits(3,2,cellLineIter);
    
    
    
    %R2(cellLineIter) = runCompartmentModel_r2(parm,concMatrix,...
    %    allTimes_hrs_mod, analysisWells, aifDuration,startDrug,endDrug);
    [Y_meas,Y_pred,Y_weight] = runCompartmentModel_CellLine_forResidCalc(parm,concMatrix,...
        allTimes_hrs_mod, analysisWells, aifDuration,startDrug,endDrug);
    
    figure(1)
    subplot(2,4,cellLineIter)
    allRes = (Y_meas-Y_pred) .* Y_weight;%./Y_meas;
    allRes(Y_meas<=0) = nan;
    plotME = nanmean(allRes,1);
    errorME = nanstd(allRes,1);
    %plot(allTimes_hrs_mod, plotME)
    %errorbar(allTimes_hrs_mod, plotME, errorME,'k-')
    %plot(allTimes_hrs_mod,plotME,'k.')
    
    
    
    dw_statistic = nansum((allRes(:,2:end)-allRes(:,1:end-1)).^2,2)./nansum(allRes.^2,2);
    min(dw_statistic)
    max(dw_statistic)
    
    hold on
    plot(-1,-1,'r.')
    plot(-1,-1,'k.')
    plot(allTimes_hrs_mod',allRes(1:8,:),'r.')
    plot(allTimes_hrs_mod',allRes(9:16,:),'k.')
    if cellLineIter ==4
        legend('6 hr treatment','12 hr treatment')
    end
    %title(CellLine)
    axis([0 24 -.03 .06])
    xlabel('Time (hours)')
    if cellLineIter ==1
        ylabel('Residuals (AUC-normalized)')
    end
    set(gca,'FontSize',14)
    
    
    %figure(2)
    subplot(2,4,cellLineIter+4)
    hold on
    t6 = allRes(1:8,:);
    t12 = allRes(9:16,:);
    [N,X] = hist(t6(~isnan(t6)), 50);
    [N2] = hist(t12(~isnan(t12)), X);
    b = bar(X,[N./sum(N);N2./sum(N2)]','grouped');
    b(1).EdgeColor = [1 0 0];
    b(1).FaceColor = [1 0 0];
    b(2).EdgeColor = [0 0 0];
    b(2).FaceColor = [0 0 0];
    b(1).BarWidth = 1;
    b(2).BarWidth = 1;
    %bar(X,N./sum(N),'r')
    %bar(X,N2./sum(N2),'k')
    %title(CellLine)
    xlabel('Residuals (AUC-normalized)')
    if cellLineIter ==4
        legend('6 hr treatment','12 hr treatment')
    end
    if cellLineIter ==1
        ylabel('Frequency (%)')
    end
    set(gca,'FontSize',14)
    
    figure(3)
    allConcs = allConcs([1:8 11:18]);
    aifDuration = aifDuration([1:8 11:18]);
    uconcs_pl = unique(allConcs);
    clrs = {'r.','k.'};
    N = cell(8,1);
    cnt = 0;
    legs = cell(4,1);
    for uci = 1:length(uconcs_pl)
        for ad = [6 12]
            cnt = cnt+1;
            cr = (allRes(allConcs == uconcs_pl(uci) & aifDuration==ad, :));
            if uci == 1
                [N{cnt},X] = hist(cr(~isnan(cr)),linspace(-.02, .04,50));
                %[N{cnt},X] = hist(cr(~isnan(cr)),50);
            else
                [N{cnt}] = hist(cr(~isnan(cr)),X);
            end
            %semilogx(uconcs_pl(uci), ,clrs{cnt})
        end
        legs{uci} = sprintf('%2.0d nM',round(uconcs_pl(uci)));
    end
    
    rawr = [];
    for iter = 1:2:length(N)
        rawr = [rawr [N{iter}./sum(N{iter})]'];
    end
    
    figure(3)
    subplot(2,4,cellLineIter)
    hold on
    clrs = distinguishable_colors(4);
    mrks = {'k-.','k:','k-','k--'};
    for iter = 1:4
        cp = plot(X,smooth(rawr(:,iter),5),mrks{iter},'LineWidth',2);
        set(cp,'Color',clrs(iter,:));
    end
    %b = bar(X,rawr,'grouped');
    axis([-.02 .04 0 .2])
    xlabel('Residuals (AUC-normalized)')
    if cellLineIter == 4
        legend(legs)
    end
    if cellLineIter == 1
        ylabel('Frequency (%)')
    end
    set(gca,'FontSize',14)
    
    rawr = [];
    for iter = 2:2:length(N)
        rawr = [rawr [N{iter}./sum(N{iter})]'];
    end
    
    subplot(2,4,cellLineIter+4)
    hold on
    clrs = distinguishable_colors(4);
    mrks = {'k-.','k:','k-','k--'};
    for iter = 1:4
        cp = plot(X,smooth(rawr(:,iter),5),mrks{iter},'LineWidth',2);
        set(cp,'Color',clrs(iter,:));
    end
    %b = bar(X,rawr,'grouped');
    axis([-.02 .04 0 .2])
    xlabel('Residuals (AUC-normalized)')
    %legend(legs)
    if cellLineIter == 1
        ylabel('Frequency (%)')
    end
    set(gca,'FontSize',14)
    
    
    
    figure(4);
    uconcs_pl = unique(allConcs);
    clrs = {'r.','k.'};
    N = cell(8,1);
    legs = cell(4,1);
    for uci = 1:length(uconcs_pl)
        subplot(4,4,cellLineIter + (uci-1)*4)
        hold on
        plot(-1,-1,'r.')
        plot(-1,-1,'k.')
        cr = (allRes(allConcs == uconcs_pl(uci) & aifDuration==6, :));
        plot(allTimes_hrs_mod',cr,'r.')
        
        cr = (allRes(allConcs == uconcs_pl(uci) & aifDuration==12, :));
        plot(allTimes_hrs_mod',cr,'k.')
        axis([0 24 -.03 .06])
    end
    
    if cellLineIter == 4
        subplot(4,4,4)
        legend('6 hr treatment','12 hr treatment','Location','Best')
    end
    
    
    
    
end
