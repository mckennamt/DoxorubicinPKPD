% full script to test fitting to simulated data with various noise levels

clear all;clc

%% generate simulated data


cellData = loadCellLineData('SUM149');

% Only create 1 dataset for simulation
cellData.ExposureTimes = 6;
cellData.times = cellData.times(1);
cellData.drugConc = cellData.drugConc(1);
cellData.Tspan = cellData.Tspan(1);
cellData.allYall = cellData.allYall(1);
cellData.allYd_tv = cellData.allYd_tv(1);
cellData.dox_pss = cellData.dox_pss(1);
cellData.uconcs = cellData.uconcs(1);

uconcs = cellData.uconcs{1};
Tspan = cellData.Tspan{1};

% define compartment parameters and kp, theta
compartmentParms = [.06 .06/.13 .05];

% [kp, theta, functionSelect]
kp = .02;
fixedParms = [.02 2e4 5];


%%



allTrue = cell(5,1);
allPred = cell(5,1);
allError = cell(5,1);
allCISize = cell(5,1);

nli = 0;
for noiseLevels = [0 .05 .1 .15 .2]
    nli = nli+1;
    holdTrue = [];
    holdPred = [];
    CISize = [];
    errorPred_percent = [];
    for numIterAtNoiseLevel = 1:2500%500
        if numIterAtNoiseLevel< 2500
            fprintf('%2d ', numIterAtNoiseLevel)
        else
            fprintf('%2d\n\n\n\n\n', numIterAtNoiseLevel)
        end
        % supply parameters for the simulations
        thetaSim = unifrnd(1.6e4, 2.4e4, 1, length(uconcs));
        kd_dum = unifrnd(0, 2*.02, 1, length(uconcs));
        r_dum = unifrnd(.001, .05, 1, length(uconcs));
        
        %specified_fitParms = [kd_dum' r_dum' thetaSim'];
        specified_fitParms = [kd_dum' thetaSim'];
        specified_fitParms(1,1) = 0;
        specified_fitParms(1,2) = 2e4;
        holdTrue = [holdTrue;specified_fitParms(2:end,:)];
        
        % run simulations
        %figure(1);clf
        for drugIter = 1:length(uconcs)
            
            fitParms = specified_fitParms(drugIter,:);
            S0 = 500 + 0*randn(1);
            
            % run simulation
            %Ypr = damageModel_analytic(S0, Tspan, fitParms, [.02 2e4 5]);
            Ypr = damageModel_analytic(S0, Tspan, fitParms, [.02 2e4 1]);
            
            % add in noise to simulations
            Ypr_sim = zeros(6,length(Ypr));
            for noisy = 1:6
                Ypr_sim(noisy,:) = Ypr + noiseLevels.*Ypr.*randn(size(Ypr));
            end
            
            Ypr_sim(Ypr_sim<=1) = 1;
            
            %figure(1);
            %hold on
            %plot(Tspan, Ypr_sim(:,:))
            %pause(.01)
            cellData.allYall{1}{drugIter} = Ypr_sim';
            
            
        end
        
        %cellData = rmfield(cellData,{'drugConc','NumNuclei','NumNuclei2',...
        %    'NucDensity','DrugTime','times'});
        
        %% fit simulated data
        
        %funToFit = 5;
        funToFit = 1;
        
        % Get fits for growth rate, theta
        % Fit to control data and pre-treatment timepoints
        % use all exposure time data sets for this
        lsqOpts = optimoptions('lsqnonlin','MaxFunEvals',10000,'Display','none',...
            'TolFun', 1e-50,'MaxIter', 5000,'TolX', 1e-5);
        prob = createOptimProblem('lsqnonlin','x0',[0 0],...
            'objective',@(x) runGrowthModel_allData(x,cellData),...
            'lb',[0 1],'ub',[1 1e6],'options',lsqOpts);
        ms = MultiStart('UseParallel',true,'Display','off');
        
        % create parameter estimates
        parmEsts{1} = linspace(.01, .05, 5);
        parmEsts{2} = linspace(5000, 1e5, 5);
        hldGrid = cell(1,numel(parmEsts));
        [hldGrid{:}] = ndgrid(parmEsts{:});
        hldGrid = cellfun(@(x) x(:), hldGrid,'UniformOutput',false);
        startPointList = cat(2,hldGrid{:});
        custpts = CustomStartPointSet(startPointList);
        [fitParms_kp,~, ~, ~,solutions] = run(ms,prob,custpts);
        
        %parms_est = [.025 10000];
        [fitParms_kp,~,res,~,~,~,J] = ...
            lsqnonlin(@(x)runGrowthModel_allData(x,cellData),fitParms_kp,...
            [], [], lsqOpts);%[0 1], [10 50000]
        fitParms_kp_ci = nlparci(fitParms_kp,res,'jacobian',J);
        kp = fitParms_kp(1); %hr^-1
        theta = fitParms_kp(2); %cell count
        
        % fit rest of model
        numFunctionsFit = length(funToFit);
        numExpTime = length(cellData.ExposureTimes);
        
        collectAICs = cell(numExpTime,numFunctionsFit);
        parameterFits = cell(numExpTime,numFunctionsFit);
        parameterCI = cell(numExpTime,numFunctionsFit);
        
        for expTimeIter = 1:numExpTime
            
            % remove control data (no drug added)
            currentData = cellData.allYall{expTimeIter};
            currentData = currentData(2:end,:);
            
            allTspan = cellData.Tspan(expTimeIter);
            allTspan = repmat(allTspan,size(currentData,1),1);
            
            for functIter = 1:numFunctionsFit
                
                % fit current model
                fixedParms = [kp theta funToFit(functIter)];
                postTreat_TP = cellData.DrugAdded_Tp + 1;
                
                % fit all concentrations at given exposure time collectively
                [allParms, allParms_ci, ~, allAICs, numParms] = ...
                    fitEachTreatmentModel_regularized(...
                    currentData, fixedParms, postTreat_TP, allTspan, 0, [],[]);
                
                holdPred = [holdPred;allParms];
                
                collectAICs{expTimeIter,functIter} = allAICs;
                parameterFits{expTimeIter,functIter} = allParms;
                parameterCI{expTimeIter,functIter} = allParms_ci;
                
                errorPred_percent = [errorPred_percent;...
                    100*abs(specified_fitParms(2:end,:) - allParms)./abs(allParms)];
                
                
                
                CISize = [CISize;...
                    [squeeze(allParms_ci(:,2,:)) - squeeze(allParms_ci(:,1,:))]'];
                
            end
        end
%         %% Display simulation fits
%         
%         figure(2);clf
%         subplot(131)
%         plot(r_dum(2:end), allParms(1:end,1),'r.')
%         hold on
%         lb = allParms(1:end,1) - squeeze(allParms_ci(1,1,1:end));
%         ub = squeeze(allParms_ci(1,2,1:end))-allParms(1:end,1);
%         errorbar(r_dum(2:end), allParms(1:end,1), ...
%             lb, ub, 'ro')
%         plot(r_dum(2:end), kd_dum(2:end),'k.')
%         title('kd')
%         
%         subplot(132)
%         plot(r_dum(2:end), allParms(1:end,2),'ro')
%         hold on
%         lb = allParms(1:end,2) - squeeze(allParms_ci(2,1,1:end));
%         ub = squeeze(allParms_ci(2,2,1:end))-allParms(1:end,2);
%         errorbar(r_dum(2:end), allParms(1:end,2), ...
%             lb, ub, 'ro')
%         plot(r_dum(2:end), r_dum(2:end),'k.')
%         title('r')
%         subplot(133)
%         plot(1:9, allParms(1:end,3),'ro')
%         hold on
%         lb = allParms(1:end,3) - squeeze(allParms_ci(3,1,1:end));
%         ub = squeeze(allParms_ci(3,2,1:end))-allParms(1:end,3);
%         errorbar(1:9, allParms(1:end,3), ...
%             lb, ub, 'ro')
%         plot(1:9, thetaSim(2:end),'k.')
%         %axis([.8*min(r2_dum) 1.1*max(r2_dum) .8*min(r2_dum) 1.1*max(r2_dum)])
%         title('theta')
%         hold off
%         
%         %% view fit timecourses
%         
%         figure(100);clf
%         for drugIter = 1:length(uconcs)-1
%             
%             
%             fitParms = allParms(drugIter,:);
%             S0 = 500 + 0*randn(1);
%             
%             %[Tpr, Ypr] = ode23(@(t,y) damageModel_Simplified_TimeMod(t,y,...
%             %    fitParms,fixedParms,[],[],[]),Tspan,S0);
%             Ypr = damageModel_analytic(S0, Tspan, fitParms, fixedParms);
%             
%             Ypr_sim = zeros(6,length(Ypr));
%             for noisy = 1:6
%                 Ypr_sim(noisy,:) = Ypr;% + .2.*Ypr.*randn(size(Ypr));
%             end
%             
%             %figure(1);clf;
%             %plot(Tpr,Ypr_sim,'k.-')
%             %pause(.1)
%             
%             Ypr_sim(Ypr_sim<=1) = 1;
%             
%             cellData.allYall_Pred{1}{drugIter} = Ypr_sim';
%             
%             hold on
%             plot(cellData.allYall_Pred{1}{drugIter},'k-')
%             plot(cellData.allYall{1}{drugIter+1},'r-')
%             
%         end
%         
%         pause(.01)
%         
        
    end
    
    
    allTrue{nli} = holdTrue;
    allPred{nli} = holdPred;
    allError{nli} = errorPred_percent;
    allCISize{nli} = CISize;
    
end
%


save('NoiseSimulations_20170327.mat')
%%

% make a matrix to define average error rates

% figure(1);clf
% 
% for nslvls = 1:nli
%     
% 
%     holdTrue = allTrue{nslvls};
%     %errorPred_percent = allError{nslvls};
%     errorPred_percent = allCISize{nslvls};
%     
% gridSize = 5;
% 
% kdRange = linspace(min(holdTrue(:,1)), max(holdTrue(:,1)), gridSize);
% rRange = linspace(min(holdTrue(:,2)), max(holdTrue(:,2)), gridSize);
% 
% errorRates_kd = zeros(gridSize-1);
% errorRates_r = zeros(gridSize-1);
% for xi = 1:(gridSize-1)
%     for yi = 1:(gridSize-1)
%         
%         ss_kp = holdTrue(:,1)>= kdRange(xi) & holdTrue(:,1) < kdRange(xi+1);
%         ss_r = holdTrue(:,2)>= rRange(yi) & holdTrue(:,2) < rRange(yi+1);
%         
%         errorRates_kd(xi,yi) = mean(abs(errorPred_percent(ss_kp&ss_r, 1)));
%         errorRates_r(xi,yi) = mean(abs(errorPred_percent(ss_kp&ss_r, 2)));
%         
%     end
% end
% 
% %errorRates_kd(isnan(errorRates_kd)) = 1;
% 
% 
% subplot(2,nli, nslvls);
% [X,Y] = meshgrid(rRange(2:end),kdRange(2:end));
% imagesc(rRange, kdRange, errorRates_kd)
% caxis([0 .05])
% xlabel('r');ylabel('kd')
% title('kd error')
% 
% subplot(2,nli, nslvls+nli);
% imagesc(rRange, kdRange, errorRates_r)
% title('r error')
% caxis([0 .05])
% xlabel('r');ylabel('kd')
% end
