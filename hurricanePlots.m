

%close all;
clear all;clc

allCellLine = {'MDAMB231','SUM149','MDAMB468'};
cellLineSelect = 2;
CellLine = allCellLine{cellLineSelect};

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
allFitParms_ci = cell(3,5);

% n = number of samples
% resnorm = sum of squares of residuals
% k = number of parameters
% AIC = n*log(resnorm/n)+(2*k)+(((2*k)*(k+1))/(n-k-1));
AICs = cell(3,5);

extc = 0;
for ExposureTime = allExposureTimes
    extc = extc+1;
    for functionSelectIter = 1:5
        
        load(genFlName(CellLine,ExposureTime, functionSelectIter))
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
        reshaped_ci = reshape(allParms_ci(:),size(allParms_ci,1)*size(allParms_ci,2),[])';
        collectTrainingParms_ci = sortrows([uc reshaped_ci(fc:end,:)],1);
        
        %%remove noisy fits
        %collectTrainingParms(collectTrainingParms(:,2)>=49,:) = [];
        
        allFitParms{extc,functionSelectIter} = collectTrainingParms;
        allFitParms_ci{extc,functionSelectIter} = collectTrainingParms_ci;
        
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

c1 = distinguishable_colors(2);


tsc = 0;
for testSet = allExposureTimes
    
    tsc = tsc+1;
    ExposureTime = testSet;
    cf = figure(ExposureTime);clf
    set(cf,'Position',[1300 100 1400 800])
    
    fsi = 0;
    for functionSelect = [1 5]
        fsi = fsi+1;
        
        load(genFlName(CellLine,ExposureTime, functionSelect))
        
        
        % AUC or CB (bound drug) or just concentration (no time considered)?
        %uc = cellfun(@(x) max(x), allYd_tv(fc:end,1));%CB
        uc = cellfun(@(x) max(x), allYd_tv(2:end,3))*ExposureTime;%AUC
        %uc = uconcs(fc:end);% just concentration
        
        for drugs = 2:length(uconcs)
            cd = uc(drugs-1);
            
            % initial cell count for forward model
            dm = mean(cellData.NumNuclei(cellData.drugConc==uconcs(drugs),5));
            
            Yd = allYd_tv{drugs,1};
            time_vector = allYd_tv{drugs,2};
            maxDrugConc = max(Yd);
            Tspan = cellData.times-cellData.times(cellData.DrugAdded_Tp);
            Tspan = Tspan(cellData.DrugAdded_Tp+1:end);
            Tspan = [Tspan Tspan(end)+24:2:(Tspan(end)+1*2)];
            
            S0 = dm;
            
            fitParms = allFitParms{tsc,fixedParms(3)};
            fitParms = fitParms(drugs-1,2:end);
            
            fitParms_Range = allFitParms_ci{tsc,fixedParms(3)};
            fitParms_Range = fitParms_Range(drugs-1,2:end);
            
            fit_fixedParms = allFixedParms{tsc,fixedParms(3)};
            fit_fixedParms(3) = fixedParms(3);
            disp([cd fixedParms(3)])
            
            
            [Tpr2, Ypr2] = ode45(@(t,y) damageModel_Simplified_TimeMod(t,y,...
                fitParms,fit_fixedParms,maxDrugConc,Yd,time_vector),Tspan,S0);
            
            figure(ExposureTime);
            subplot(3,3,drugs-1)
            
            hold on
            
            
            
            fp_len = length(fitParms);
            fpr_len = length(fitParms_Range);
            
            holdParm = cell(1,fp_len);
            for genSamp = 1:fp_len
                holdParm{genSamp} = fitParms_Range(genSamp:fp_len:end);
            end
            allCombinations = cell(1, numel(holdParm));
            [allCombinations{:}] = ndgrid(holdParm{:});
            allCombinations = cell2mat( cellfun(@(v)v(:), allCombinations, 'UniformOutput',false));
            
            tpr_hurricane = cell(size(allCombinations,1),1);
            ypr_hurricane = cell(size(allCombinations,1),1);
            
            if functionSelect==1
                allCombinations(allCombinations(:,2)<0,2) = 1e4;
            end
            
            if functionSelect==5
                allCombinations(allCombinations(:,2)<0,2) = 0;
                allCombinations(allCombinations(:,3)<1e4,3) = 1e4;
            end
            
            if functionSelect==4
                allCombinations(allCombinations(:,2)<0,2) = 0;
                allCombinations(allCombinations(:,1)<0,1) = 0;
                allCombinations(allCombinations(:,4)<0,4) = 1e4;
            end
            
            Tspan_mod = linspace(Tspan(1),Tspan(end),100);
            for comboIter = 1:size(allCombinations,1)
                [tpr_hurricane{comboIter}, ypr_hurricane{comboIter}] = ode23s(@(t,y) damageModel_Simplified_TimeMod(t,y,...
                    allCombinations(comboIter,:),fit_fixedParms,maxDrugConc,Yd,time_vector),Tspan_mod,S0);
                %plotMe = ypr_hurricane{comboIter};
                %currplot = plot(tpr_hurricane{comboIter},plotMe(:,1),'k-','LineWidth',2);
                %set(currplot,'Color',c1(fsi,:))
            end
            
            
            allProjections = [ypr_hurricane{:}];
            theShade = zeros(length(Tspan_mod)*2,2);
            theShade(:,1) = [Tspan_mod'; Tspan_mod(end:-1:1)'];
            for makeShade = 1:length(Tspan_mod)
                theShade(makeShade,2) = min(allProjections(makeShade,:));
            end
            shadeIndex = length(Tspan_mod);
            for makeShade = length(Tspan_mod):-1:1
                shadeIndex = shadeIndex+1;
                theShade(shadeIndex,2) = max(allProjections(makeShade,:));
            end
            
            
            theAlphaVal = .1;
            patch(theShade(:,1),theShade(:,2),c1(fsi,:),'FaceAlpha',theAlphaVal,...
                'EdgeColor',c1(fsi,:),'EdgeAlpha',.2,'LineStyle','-')
            
            
            currplot = plot(Tpr2,Ypr2(:,1),'k-','LineWidth',2);
            set(currplot,'Color',c1(fsi,:))
            
            %legend('Prediction','Best Fit','Location','Best')
            hold off
            
        end
        
        figure(ExposureTime)
        for drugs = 2:length(uconcs)%cd%cd-1:cd+1%
            subplot(3,3,drugs-1)
            hold on
            mean_normcellData = (cellData.NumNuclei(...
                cellData.drugConc==uconcs(drugs),:));
            %for a = 1:size(mean_normcellData,1)
            %    plot(cellData.times(:)-cellData.times(cellData.DrugAdded_Tp),mean_normcellData(a,:),'k.',...
            %        'MarkerSize',15,'LineWidth',2)
            %end
            
            daStd = std(mean_normcellData);
            daMean = mean(mean_normcellData);
            nSamp = size(mean_normcellData,1);
            
            bnds = tinv([0.025, 0.975],nSamp-1)' * daStd/sqrt(nSamp);
            ub = bnds(2,:);
            lb = -bnds(1,:);
            
            %[mean_normcellData,~,meanCI] = normfit(cellData.NumNuclei(...
            %    cellData.drugConc==uconcs(drugs),:));
            %ub = meanCI(2,:)-mean_normcellData;
            %lb = mean_normcellData - meanCI(1,:);
            
            %mean_normcellData2 = mean(cellData.NumNuclei(...
            %    cellData.drugConc==uconcs(drugs),:));
            %ub2 = max(cellData.NumNuclei(...
            %    cellData.drugConc==uconcs(drugs),:)) - mean_normcellData;
            %lb2 = mean_normcellData - min(cellData.NumNuclei(...
            %    cellData.drugConc==uconcs(drugs),:));
            
            errorbar(cellData.times(:)-cellData.times(cellData.DrugAdded_Tp),daMean,...
                lb,ub,...
                'ks','MarkerSize',2,'LineWidth',2)
            
            %errorbar(cellData.times(:)-cellData.times(cellData.DrugAdded_Tp),mean_normcellData2(1,:),...
            %    lb2,ub2,...
            %    'cs','MarkerSize',2,'LineWidth',2)
            
            daTitle = sprintf('%3.0f nM',uconcs(drugs));
            title(daTitle)
            
            %set(gca,'XLim',[-120 (Tspan(end)+1*24)])
            set(gca,'XLim',[-120 600])
            set(gca,'YLim',[0 2.5e4])
        end
        %pause(.1)
        %xlabel('Time (hours)')
        %ylabel('Population Size')
        
    end
end
