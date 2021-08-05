

%close all;clear all;clc

% Leave one out analysis to train model and make predictions
% Doesn't predict new concentration at given exposure time, only similar
% AUC's


% which function to use?
% 1 = single parameter
% 2 = 'linear' model
% 3 = 'exponential' model
% 4 = 'biexponential' model
% 5 = 'single biexponential' model

%% load all parameter fits
allExposureTimes = [6 12 24];
allFixedParms = cell(3,5);
allFitParms = cell(3,5);

% n = number of samples
% resnorm = sum of squares of residuals
% k = number of parameters
% AIC = n*log(resnorm/n)+(2*k)+(((2*k)*(k+1))/(n-k-1));
AICs = cell(3,5);

extc = 0;
for ExposureTime = 6
    extc = extc+1;
    for functionSelectIter = 5
        
        allFixedParms{extc,functionSelectIter} = fixedParms;
        
        % first concentration is untreated, useless parameter
        fc = 2;
        
        % AUC or CB (bound drug) or just concentration (no time considered)?
        %uc = cellfun(@(x) max(x), allYd_tv(fc:end,1));%CB
        uc = cellfun(@(x) max(x), allYd_tv(fc:end,3))*ExposureTime;%AUC
        %         %uc = uconcs(fc:end);% just concentration
        
        % first column is concentration data, subsequent columns parameter 1, 2,...
        % sort by concentration
        collectTrainingParms = sortrows([uc allParms(fc:end,:)],1);
        
        %%remove noisy fits
        %collectTrainingParms(collectTrainingParms(:,2)>=49,:) = [];
        
        allFitParms{extc,functionSelectIter} = collectTrainingParms;
        
        collectAIC = zeros(length(uconcs)-1,2);
        k = numParms+1;%estimating sigma^2 too!
        
        for drugIter = 2:length(uconcs)
            
            gParms = allParms(drugIter,:);
            Yd = allYd_tv{drugIter,1};
            time_vector = allYd_tv{drugIter,2};
            maxDrugConc = max(Yd);
            
            cYall = allYall{drugIter};
            cYall = cYall(postTreat_TP:end,:);
            Tspan_modified = Tspan(postTreat_TP:end);
            
        end
        
        
    end
end


%% Run prediction forward, compare to ground truth

% add in chi-square goodness of fit calculation (compare to fit)
% xi_sq= sum((observed-expected).^2/expected);

% compAIC = [AUC value, AIC for predicted model, AIC for best fit,
% delta(AIC)]

c1 = distinguishable_colors(5);


tsc = 0;
for testSet = 6
    
    tsc = tsc+1;
    ExposureTime = testSet;
    cf = figure(ExposureTime);clf
    set(cf,'Position',[100 100 900 500])
    
    for functionSelect = 5
        
        
        % AUC or CB (bound drug) or just concentration (no time considered)?
        %uc = cellfun(@(x) max(x), allYd_tv(fc:end,1));%CB
        uc = cellfun(@(x) max(x), allYd_tv(2:end,3))*ExposureTime;%AUC
        %uc = uconcs(fc:end);% just concentration
        
        for drugs = 2:length(uconcs)
            cd = uc(drugs-1);
            
            % initial cell count for forward model
            dm = mean(cellData.NumNuclei(cellData.drugConc==uconcs(drugs),cellData.DrugAdded_Tp+1));
            
            Yd = allYd_tv{drugs,1};
            time_vector = allYd_tv{drugs,2};
            maxDrugConc = max(Yd);
            Tspan = cellData.times-cellData.times(cellData.DrugAdded_Tp);
            Tspan = Tspan(cellData.DrugAdded_Tp+1:end);
            Tspan = [Tspan Tspan(end)+24:24:(Tspan(end)+1*24)];
            
            S0 = dm;
            
            fitParms = allFitParms{tsc,fixedParms(3)};
            fitParms = fitParms(drugs-1,2:end);
            
            fit_fixedParms = allFixedParms{tsc,fixedParms(3)};
            fit_fixedParms(3) = fixedParms(3);
            disp([cd fixedParms(3)])
            
            
            [Tpr2, Ypr2] = ode45(@(t,y) damageModel_Simplified_TimeMod(t,y,...
                fitParms,fit_fixedParms,maxDrugConc,Yd,time_vector),Tspan,S0);
            
            
            cYall = allYall{drugs};
            cYall = cYall(postTreat_TP:end,:);
            Tspan_modified = cellData.times-cellData.times(cellData.DrugAdded_Tp);
            Tspan_modified = Tspan_modified(postTreat_TP:end);
            
            [resNorm_bf,n_bf,percentDiff_bf] = runDamageModel_AICCalc(fitParms,cYall,Tspan_modified,...
                fit_fixedParms,maxDrugConc,Yd,time_vector);
            
            
            figure(ExposureTime);
            subplot(3,3,drugs-1)
            hold on
            currplot = plot(Tpr2,Ypr2(:,1),'k-','LineWidth',2);
            set(currplot,'Color',c1(functionSelect,:))
            %legend('Prediction','Best Fit','Location','Best')
            hold off
            
        end
        
        figure(ExposureTime)
        for drugs = 2:length(uconcs)%cd%cd-1:cd+1%
            subplot(3,3,drugs-1)
            hold on
            mean_normcellData = cellData.NumNuclei(...
                cellData.drugConc==uconcs(drugs),:);
            for a = 1:size(mean_normcellData,1)
                plot(cellData.times(:)-cellData.times(cellData.DrugAdded_Tp),mean_normcellData(a,:),'kp',...
                    'MarkerSize',2,'LineWidth',2)
            end
            mean_normcellData = mean(cellData.NumNuclei(...
                cellData.drugConc==uconcs(drugs),:));
            daTitle = sprintf('%3.0f nM',uconcs(drugs));
            title(daTitle)
            
            set(gca,'XLim',[-120 (Tspan(end)+1*24)])
        end
        
        %xlabel('Time (hours)')
        %ylabel('Population Size')
        
    end
end
%%

figure(1212);clf

% first concentration is untreated, useless parameter
fc = 2;

% AUC or CB (bound drug) or just concentration (no time considered)?
%uc = cellfun(@(x) max(x), allYd_tv(fc:end,1));%CB
uc = cellfun(@(x) max(x), allYd_tv(fc:end,3))*ExposureTime;%AUC
%uc = uconcs(fc:end);% just concentration


for iter = 1:numParms
    cp = allParms(2:end,iter);
    
    lb = squeeze(allParms_ci(iter,1,2:end));
    ub = squeeze(allParms_ci(iter,2,2:end));
    
    subplot(1,numParms,iter)
    semilogx(uc,cp,'r.','MarkerSize',1,'LineWidth',.1)
    hold on
    errorbar(uc,cp,abs(cp-lb),abs(ub-cp),'r.','LineWidth',.2)
    
	semilogx(uc, specified_fitParms(2:end,iter),'ko')

    
end
