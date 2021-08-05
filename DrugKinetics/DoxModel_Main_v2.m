
% Dynamic dox modeling
close all;clear all;clc

CellLine = 'MDAMB453';
SelectCompartmentModelExperiment

load(fullfile(pwd, ['SegmentationMasks_' CellLine '_complete.mat']));
%load(fullfile(pwd, ['SegmentationMasks_' CellLine '_20160914_complete.mat']));

% all the below information should have been read in during the cell
% segmenation process
% Do it anyway!

% Read in plate map file to get concentrations and wells of interest
isDynamic = true;
[plateMap,DrugTime,cellType,folderList] = importPlateMap(PlateMapFile,isDynamic);

% Wells of interest...
% Define control wells
drugOnlyWells = cell2mat(cellfun(@(x) strcmp(x,'drugOnly'),cellType,'UniformOutput',false));
drugOnlyFolders = cell(sum(drugOnlyWells),1);
iter = 0;
for cwells = 1:length(drugOnlyWells)
    if drugOnlyWells(cwells)
        iter = iter+1;
        drugOnlyFolders{iter} = ['Well ' folderList{cwells}];
    end
end
drugOnlyConcs = plateMap{3}(drugOnlyWells);

% Define aif wells
aifWells = cell2mat(cellfun(@(x) strcmp(x,'aif'),cellType,'UniformOutput',false));
aifFolders = cell(sum(aifWells),1);
iter = 0;
for cwells = 1:length(aifWells)
    if aifWells(cwells)
        iter = iter+1;
        aifFolders{iter} = ['Well ' folderList{cwells}];
    end
end
aifConcs = plateMap{3}(aifWells);
aifDuration = plateMap{4}(aifWells);

% Wells with cells
pbsWells = cell2mat(cellfun(@(x) strcmp(x,'null'),cellType,'UniformOutput',false));
cellWells = ~(pbsWells | drugOnlyWells | aifWells);
cellFolders = cell(sum(cellWells),1);
iter = 0;
for cwells = 1:length(cellWells)
    if cellWells(cwells)
        iter = iter+1;
        cellFolders{iter} = ['Well ' folderList{cwells}];
    end
end
cellDrugConcs = plateMap{3}(cellWells);
cellDrugDuration = plateMap{4}(cellWells);

% All folders of interest
foi = [drugOnlyFolders; aifFolders; cellFolders];
allConcs = cellDrugConcs;

analysisWells = [1:8 11:18];

%% Get time vectors and number of images

% Should be loaded with cell segmentation data
c_foi = foi{1};
[allTimes_hrs, numImage_perFolder] = createTimeVector(DataFolders,c_foi);
allTimes_hrs = allTimes_hrs - allTimes_hrs(startDrug+1);

%% Read and process images, single exposure time

% fluorIntenMatrix = number of wells x number of images x 4
% 4: mean cell value, mean in AIF, std cell, std AIF
%fluorIntenMatrix = zeros(1,sum(numImage_perFolder),4);
%fluorIntenMatrix = zeros(8,50,4);

% processing
% subtract off corresponding aif well from cell well
% Signal is summation of 2 compartments (cytoplasm, nuclear)

maxInten = 1;%1000
davec = 0:.01:1; %used for histogram correction


for wellIter = analysisWells
    disp(wellIter)
    
    for dfIter = 1:length(DataFolders)
        
        controlFolder_blank = fullfile(DataFolders{dfIter},aifFolders{9});
        controlFolder_blank2 = fullfile(DataFolders{dfIter},aifFolders{10});
        
        controlFolder_AIF = fullfile(DataFolders{dfIter},aifFolders{wellIter});
        cellFolder = fullfile(DataFolders{dfIter},cellFolders{wellIter});
        drugOnlyFolder = fullfile(DataFolders{dfIter},drugOnlyFolders{wellIter});
        
        allOuttie =  zeros(numImage_perFolder(dfIter),4);
        fewImages = numImage_perFolder(dfIter)-1;
        fewImages = fewImages-1;
        
        for allIm = 0:fewImages
            
            % Need to re-write processDynamicFluorescentImages to do the
            % histogram normalization technique
            % correct images before processing
            
            flName = sprintf('Doxorubicin - n%06d.tif',allIm);
            
            aif_blank_flName = fullfile(controlFolder_blank,flName);
            aif_blank2_flName = fullfile(controlFolder_blank,flName);
            aif_flName = fullfile(controlFolder_AIF,flName);
            cell_flName = fullfile(cellFolder,flName);
            
            % wells with drug added
            Im_cell = double(imread(cell_flName));
            Im_aif = double(imread(aif_flName));
            
            % average background signal (no drug)
            Im_blank = double(imread(aif_blank_flName));
            Im_blank2 = double(imread(aif_blank2_flName));
            Im_blank = (Im_blank+Im_blank2)./2;
            
            % get cell mask from cell segmentation method
            cCellMask = CellMasks{wellIter,allIm+1};
            
            %cvals = Im_cell(cCellMask);mcv = mean(cvals);
            %avals = Im_aif(cCellMask);
            %dDiff = cvals-avals;
            %mean(cvals(dDiff>=0))
            %mean(Im_cell(cCellMask))
            %mean(cvals(cvals>mcv))
            
            c_mean = Im_cell(cCellMask);
            srted_mean = sort(c_mean);
            ioi = ceil(.05*length(srted_mean));
            ioi2 = ceil(.95*length(srted_mean));
            clip_mean = mean(srted_mean(ioi:ioi2));
            %clip_mean = mean(c_mean);
            
            %figure(87);clf;
            %subplot(121)
            %imagesc(Im_cell - Im_aif);caxis([-10 100]);
            %colormap gray
            %hold on
            %lab = label2rgb(cCellMask,'autumn');
            %h = imshow(lab);
            %set(h,'AlphaData',.1)
            %subplot(122)
            %imDiff = Im_cell - Im_aif;
            %%imDiff(imDiff<-20) = 0;
            %%imDiff(imDiff>100) = 0;
            %[n,x] = hist(imDiff(cCellMask),-20:10:100);
            %plot(clip_mean, linspace(0,1,1000),'r-','LineWidth',5)
            %hold on
            %bar(x,n./max(n),'g')
            %axis([-20 100 0 1.1])
            %pause
            
            %allOuttie(allIm+1,:) = [mean(Im_cell(cCellMask)) mean(Im_cell(~cCellMask)),...
            %    mean(Im_aif(:)) mean(Im_blank(:))];
            allOuttie(allIm+1,:) = [clip_mean mean(Im_cell(~cCellMask)),...
                mean(Im_aif(:)) mean(Im_blank(:))];
            
            
            
        end
        
        startInsert = sum(numImage_perFolder(1:dfIter-1)) + 1;
        endInsert = sum(numImage_perFolder(1:dfIter));
        fluorIntenMatrix(wellIter, startInsert:endInsert, :) = allOuttie;
        
    end
    
end

concMatrix = fluorIntenMatrix;
if strcmp(CellLine,'MDAMB453')
    concMatrix(:,54,:) = [];
    fluorIntenMatrix(:,54,:) = [];
    allTimes_hrs(54) = [];
end


%% for each timepoint, make concentration look up table

concFluorInten_Map = zeros(length(allTimes_hrs),2);

for dfIter = 1:length(DataFolders)
    fewImages = numImage_perFolder(dfIter)-1;
    
    for allIm = 0:fewImages
        fluorInten = zeros(size(drugOnlyConcs));
        stdInten = zeros(size(drugOnlyConcs));
        collectIt = zeros(1024*1344, length(drugOnlyFolders));
        for wellIter = 1:length(drugOnlyFolders)
            
            drugOnlyFolder = fullfile(DataFolders{dfIter},drugOnlyFolders{wellIter});
            
            flName = sprintf('Doxorubicin - n%06d.tif',allIm);
            drugOnly_flName = fullfile(drugOnlyFolder,flName);
            Im_drug = double(imread(drugOnly_flName));
            
            collectIt(:,wellIter) = Im_drug(:);
            
            fluorInten(wellIter) = mean(Im_drug(:));
            stdInten(wellIter) = std(Im_drug(:));
        end
        
        uConc = unique(drugOnlyConcs);
        fluorInten = zeros(size(uConc));
        stdInten = zeros(size(uConc));
        for uciter = 1:length(uConc)
            hld = collectIt(:, drugOnlyConcs == uConc(uciter));
            fluorInten(uciter) = mean(hld(:));
            stdInten(uciter) = std(hld(:));
            
        end
        
        % fit linear model to establish intensity vs concentration map at
        % each timepoint
        linFit = polyfit(uConc,fluorInten,1);
        concFluorInten_Map(allIm+1,:) = linFit;
        
        lb = -1.96*stdInten;
        ub = 1.96*stdInten;
        
        cf = figure(1);clf;
        %plot(uConc,fluorInten,'r.')
        ss = [1 4:length(uConc)];
        hold on
        cMod = linspace(0,2800);
        fMod = linFit(1)*cMod + linFit(2);
        plot(cMod,fMod,'k-','LineWidth',1.5)
        errorbar(uConc(ss), fluorInten(ss), lb(ss), ub(ss),'ro')
        axis([-100 3000 200 500])
        %set(cf,'Position',[100 100 900 250])
        set(gca,'FontSize',18)
        
        %pause
    end
    
end
%
% % convert to concentration prior to fitting model parameters
% % fluorIntenMatrix_mod = fluorIntenMatrix(:,1:fewImages+1,:);
% % %fluorIntenMatrix_mod(:, :, 2) = fluorIntenMatrix_mod(:,:, 3);
% % concMatrix = convertToConcentration(fluorIntenMatrix_mod, allConcs, allTimes_hrs);
% %
%
% concMatrix = fluorIntenMatrix;
%% fit compartment model

allTimes_hrs_mod = allTimes_hrs + abs(min(allTimes_hrs));
%options=optimset('TolFun',1e-6,'TolX',1e-6,'MaxIter',1e50,'Display','off','PrecondBandWidth',0);

estim =     [.01 .02 10];
lowBound = [0 0 0];
highBound = [1 1 100];

% 3 compartment 3 parameter
% kpars =   [k1 k2 k3]
%estim =     [.01 .02 10 .05];
%estim =     [2 0.05 .05];
%lowBound =  [0.001  0.001 .001];
%highBound = [10  5 5];
% create parameter estimates
parmEsts{1} = linspace(.01, .2, 5);
parmEsts{2} = linspace(.01, .2, 5);
parmEsts{3} = linspace(.1, 20, 5);
hldGrid = cell(1,numel(parmEsts));
[hldGrid{:}] = ndgrid(parmEsts{:});
hldGrid = cellfun(@(x) x(:), hldGrid,'UniformOutput',false);
startPointList = cat(2,hldGrid{:});
custpts = CustomStartPointSet(startPointList);



%[parm,Rn,Res,whyExit,~,~,J] = lsqnonlin(@(kpars) runCompartmentModel_CellLine(kpars,concMatrix,...
%    allTimes_hrs_mod, analysisWells, aifDuration,startDrug, endDrug),...
%    estim, lowBound, highBound, lsqOpts);

lsqOpts = optimoptions('lsqnonlin','MaxFunEvals',10000,'Display','none',...
    'TolFun', 1e-50,'MaxIter', 5000,'TolX', 1e-5);
prob = createOptimProblem('lsqnonlin','x0',startPointList(1,:),...
    'objective',@(kpars) runCompartmentModel_CellLine(kpars,concMatrix,...
    allTimes_hrs_mod, wellIter, aifDuration,startDrug, endDrug),...
    'lb',lowBound,'ub',highBound,'options',lsqOpts);
ms = MultiStart('UseParallel',true,'Display','off');
[parm,~, ~, ~,ms_solutions] = run(ms,prob,custpts);

% run at the solution to calculate the residual and jacobian to calculate
% the confidence intervals for each parameter
[parm,Rn,Res,whyExit,~,~,J] = lsqnonlin(@(kpars) runCompartmentModel_CellLine(kpars,concMatrix,...
    allTimes_hrs_mod, wellIter, aifDuration,startDrug, endDrug),...
    parm, lowBound, highBound, lsqOpts);

%[parm,Rn,Res,whyExit,~,~,J] = lsqnonlin(@(kpars) runCompartmentModel_CellLine(kpars,concMatrix,...
%    allTimes_hrs_mod, [6], aifDuration,startDrug, endDrug),...
%    estim, lowBound, highBound, options);

dof= (length(Res)-length(parm));
v = Rn/dof;
% the parameter covariance matrix
C = v*inv(J'*J);

% extract the standard errors for each fitted parameter
se_g = sqrt(diag(C));

% compute the 95% confidence intervals
tst = tinv(0.975,dof);
%disp(['95% confidence intervals of estimates of parameters'])
z_ci = [parm'-tst*se_g,parm',parm'+tst*se_g]

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% old way of fitting
% % 3 compartment 3 parameter
%     % kpars =   [k1 k2 k3]
%     estim =     [2 0.05 .05];
%     lowBound =  [0.001  0.001 .001];
%     highBound = [10  5 5];
%     [parm,Rn,Res,whyExit,~,~,J] = lsqcurvefit(@(kpars,cp) modelCt3C3P(kpars,cp,allTimes_hrs_mod),...
%         estim, inputFcn, cellUptakeSignal,...
%         lowBound, highBound, options);
%     ks = parm;
%     ks(2) = parm(1)/parm(2);
%
%     dof= (length(Res)-length(parm));
%     v = Rn/dof;
%     % the parameter covariance matrix
%     C = v*inv(J'*J);
%
%     % extract the standard errors for each fitted parameter
%     se_g = sqrt(diag(C));
%
%     % compute the 95% confidence intervals
%     tst = tinv(0.975,dof);
%     %disp(['95% confidence intervals of estimates of parameters'])
%     z_ci(:,:,cnt) = [parm'-tst*se_g,parm',parm'+tst*se_g];
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% View compartment model fit for each data well

%estim =     [.1 .1 10 .05];
%lowBound = [0 0 0 0];highBound = [1 1 Inf 1];
estim =     [.01 .02 10];
lowBound = [0 0 0];
highBound = [1 1 100];


% create parameter estimates
parmEsts{1} = linspace(.01, .2, 5);
parmEsts{2} = linspace(.01, .2, 5);
parmEsts{3} = linspace(.1, 20, 5);
hldGrid = cell(1,numel(parmEsts));
[hldGrid{:}] = ndgrid(parmEsts{:});
hldGrid = cellfun(@(x) x(:), hldGrid,'UniformOutput',false);
startPointList = cat(2,hldGrid{:});
custpts = CustomStartPointSet(startPointList);


for wellIter = analysisWells
    
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
    
    
    lsqOpts = optimoptions('lsqnonlin','MaxFunEvals',10000,'Display','none',...
        'TolFun', 1e-50,'MaxIter', 5000,'TolX', 1e-5);
    prob = createOptimProblem('lsqnonlin','x0',startPointList(1,:),...
        'objective',@(kpars) runCompartmentModel_CellLine(kpars,concMatrix,...
        allTimes_hrs_mod, wellIter, aifDuration,startDrug, endDrug),...
        'lb',lowBound,'ub',highBound,'options',lsqOpts);
    ms = MultiStart('UseParallel',true,'Display','off');
    [parm,~, ~, ~,ms_solutions] = run(ms,prob,custpts);
    
    % run at the solution to calculate the residual and jacobian to calculate
    % the confidence intervals for each parameter
    [parm,Rn,Res,whyExit,~,~,J] = lsqnonlin(@(kpars) runCompartmentModel_CellLine(kpars,concMatrix,...
        allTimes_hrs_mod, wellIter, aifDuration,startDrug, endDrug),...
        parm, lowBound, highBound, lsqOpts);
    dof= (length(Res)-length(parm));
    v = Rn/dof;
    C = v*inv(J'*J);
    se_g = sqrt(diag(C));
    tst = tinv(0.975,dof);
    z_ci = [parm'-tst*se_g,parm',parm'+tst*se_g]
    
    
    allTimes_upreg = linspace(0,24,100);
    moreTime = [allTimes_hrs_mod(end)];%+1:1:50]';
    allTimes_hrs_mod = [allTimes_hrs_mod];%; moreTime];
    inputFcn_mod = [inputFcn];% 0.*ones(1,length(moreTime))];
    %inputFcn_mod = zeros(size(inputFcn_mod));
    %inputFcn_mod(startDrug:c_endDrug+1) = max(inputFcn(:));
    upregInput = interp1(allTimes_hrs_mod,inputFcn_mod,allTimes_upreg);
    
    fwdModEval = modelCt3C3P(parm,upregInput,allTimes_upreg);
    
    %[~,fwdModEval] = ode45(@(t,y) doxCarrierModel(t,y,parm,upregInput, allTimes_upreg),...
    %    allTimes_upreg, [cellUptakeSignal(2) cellUptakeSignal(2)]);
    %fwdModEval = sum(fwdModEval,2);
    
    %[~,fwdModEval] = ode45(@(t,y) twoPeakDoxModel(t,y,parm,upregInput, allTimes_upreg),...
    %    allTimes_upreg, [cellUptakeSignal(2)]);
    
    figure(997);clf
    hold on
    plot(allTimes_hrs, (cellUptakeSignal), 'ko')
    plot(allTimes_hrs, inputFcn, 'bs')
    plot(allTimes_upreg - allTimes_hrs_mod(startDrug+1), fwdModEval, 'r-')
    plot(allTimes_upreg - allTimes_hrs_mod(startDrug+1), upregInput,'g-')
    legend('Cell uptake','Input Function','Model Fit','Location','Best')
    xlabel('Time (hrs)')
    ylabel('Concentration (nM)')
    hold off
    
end

save([CellLine 'CompartmentParameters.mat'],'concMatrix','z_ci',...
    'allTimes_hrs','allTimes_hrs_mod','analysisWells',...
    'startDrug','endDrug','aifDuration','drugOnlyConcs',...
    'allConcs','concFluorInten_Map')


