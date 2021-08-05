

cf = figure(997);clf
hold on

clrs = {'r.-','bx-'};

theLoadFls = {'MDAMB468CompartmentParameters.mat','SUM149CompartmentParameters.mat'};

for treatConds = 1:2
    for dubs = 2%:2
        load(theLoadFls{dubs})
        
        
        
        
        allTimes_hrs_mod = allTimes_hrs + abs(min(allTimes_hrs));
        if dubs == 1
            parentTimes = allTimes_hrs_mod;
        end
        
        cnt = 0;
        
        wellIter = 1;
        cnt = cnt+1;
        drugOnlyConcs(wellIter)
        options=optimset('TolFun',1e-6,'TolX',1e-6,'MaxIter',1e50,'Display','off','PrecondBandWidth',0);
        
        % Add in weighting vector?
        weight_v = ones(1,size(concMatrix,2));
        
        %concMatrix(aaa,1:28,1) = 0;
        %concMatrix(aaa,1:28,2) = 0;
        
        daif = squeeze(concMatrix(wellIter,1:end,3));
        
        inputFcn = squeeze(concMatrix(wellIter,1:end,2));
        
        blankFcn = squeeze(concMatrix(wellIter,1:end,4));
        
        if aifDuration(wellIter) == 6
            c_endDrug = endDrug(1);
        elseif aifDuration(wellIter) == 12
            c_endDrug = endDrug(2);
        end
        
        cellUptakeSignal = squeeze(concMatrix(wellIter,1:end,1));
        cellUptakeSignal(1:startDrug) = cellUptakeSignal(1:startDrug) - blankFcn(1:startDrug);
        cellUptakeSignal(c_endDrug+2:end) = cellUptakeSignal(c_endDrug+2:end) - blankFcn(c_endDrug+2:end);
        cellUptakeSignal(startDrug+1:c_endDrug+1) = cellUptakeSignal(startDrug+1:c_endDrug+1) - inputFcn(startDrug+1:c_endDrug+1);
        
        inputFcn = inputFcn-blankFcn;
        if dubs ==1
            parentInputFcn = inputFcn;
        end
        
        % 3 compartment 3 parameter
        % kpars =   [k1 k2 k3]
        estim =     [2 0.05 .05];
        lowBound =  [0.001  0.001 .001];
        highBound = [10  5 5];
        [parm,Rn,Res,whyExit,~,~,J] = lsqcurvefit(@(kpars,cp) modelCt3C3P(kpars,cp,allTimes_hrs_mod),...
            estim, inputFcn, cellUptakeSignal,...
            lowBound, highBound, options);
        ks = parm;
        ks(2) = parm(1)/parm(2);
        
        if treatConds == 1
            allTimes_upreg = linspace(0,50,100);
            upregInput = interp1(parentTimes,parentInputFcn,allTimes_upreg);
            
            upregInput = zeros(size(allTimes_upreg));
            upregInput(allTimes_upreg>3 & allTimes_upreg<9) = 100;
        else
            allTimes_upreg = linspace(0,50,100);
            upregInput = interp1(parentTimes,parentInputFcn,allTimes_upreg);
            
            upregInput = zeros(size(allTimes_upreg));
            upregInput(allTimes_upreg>3 & allTimes_upreg<6) = 200;
        end
        
        fwdModEval = modelCt3C3P(parm,upregInput,allTimes_upreg);
        hold on
        %if dubs == 1
        %    plotyy(allTimes_upreg(1:1:end) - parentTimes(5), (fwdModEval(1:1:end)),allTimes_upreg - parentTimes(5), (upregInput))
        %end
        
        p = plot(allTimes_upreg(1:1:end) - parentTimes(5), (fwdModEval(1:1:end)), clrs{dubs});
        
        xlabel('Time (hrs)')
        ylabel('Concentration (nM)')
        hold off
        
    end
    fwdModEval(end)
    hold on
    pseudoT = linspace(-5, 30, 200);
    plot(pseudoT, (ones(size(pseudoT)).*fwdModEval(end)), 'b-')
end
    %legend('Treatment Curve','MDA-MB-468','SUM-149PT','Location','Best')
    %legend('Treatment Curve','MDA-MB-468','Dox_{PSS}','Location','Best')
    %axis([-2 30 0 5])
    
    %%
    %set(cf,'Position',[1200 300 500 200])