

close all;clear all;clc

% Traditional Dose response

cf = figure(1);clf;
cColors = {'k','r','b'};
cPts = {'.','o','x'}
allExpTimes = [6 12 24];

semilogx(1:2,5:6,'k.-')
hold on
semilogx(1:2,5:6,'ro-')
semilogx(1:2,5:6,'bx-')


for cellLineIter = [3]
    allCellLine = {'MDAMB231','SUM149','MDAMB468'};
    CellLine = allCellLine{cellLineIter};
    for ExposureTimeIter = 1:3
        ExposureTime = allExpTimes(ExposureTimeIter);
        
        functionSelectIter = 1;
        
        if functionSelectIter == 1
            loadFile = sprintf('Fits_%s_SingleParm_%1dhr.mat',CellLine,ExposureTime);
        end
        
        load(loadFile)
        
        
        uconcs = unique(cellData.drugConc);
        
        evalTime = cellData.DrugAdded_Tp+1+3;
        numCells_tp8 = zeros(size(uconcs));
        
        for drugs = 1:length(uconcs)%cd-1:cd+1%
            
            %dm = mean(cellData.NumNuclei(cellData.drugConc==uconcs(drugs),1));
            %dm = mean(cellData.NumNuclei(cellData.drugConc==uconcs(drugs),4));
            numCells_tp8(drugs) = mean(cellData.NumNuclei(cellData.drugConc==uconcs(drugs),evalTime) - ...
                cellData.NumNuclei(cellData.drugConc==uconcs(drugs),cellData.DrugAdded_Tp+1));
            %dm = max(cellData.NumNuclei(cellData.drugConc==uconcs(drugs),5));
            
        end
        
        normedCount = numCells_tp8./numCells_tp8(1);
        
        
        semilogx(uconcs,normedCount,[cColors{ExposureTimeIter} cPts{ExposureTimeIter}]);
        
        
        cpFunction = @(x, drugConc) (x(1) + ((x(2)-x(1))./(1+(drugConc./x(3)).^x(4))));
        cpFunction_Fit = @(x, drugConc, Parms) (cpFunction(x,drugConc) - Parms);
        parmEst = [.01 -.05 1e3 4];
        
        [outputModelParms,~,res,~,~,~,J] = lsqnonlin(@(x) ...
            cpFunction_Fit(x,uconcs,normedCount), parmEst,[],[],lsqOpts);
        optParm_ci = nlparci(outputModelParms,res,'jacobian',J);
        
        outputModelParms(3) - optParm_ci(3,1)
        
        fprintf('%s EC50:%2.1f (%2.1f,%2.1f)\n',CellLine,outputModelParms(3),optParm_ci(3,1),optParm_ci(3,2))
        
        hold on
        modConc = logspace(0,4,100);
        fitModel = cpFunction(outputModelParms,modConc);
        semilogx(modConc,fitModel,[cColors{ExposureTimeIter} '-'])
        
        
    end
end

axis([5 3000 -.5 1.5])

legend('6 hour','12 hour','24 hour')
set(cf,'Position',[200 200 600 350])
set(gca,'FontSize',14)