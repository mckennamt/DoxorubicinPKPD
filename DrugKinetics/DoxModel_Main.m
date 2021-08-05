
% Dynamic dox modeling
close all;clear all;clc

% Define experiment folder
%ExperimentFolder = '/media/matthew/DataDisk/DrugKinetics_Data/20150727_DoxCompartmentModel_Take1';
%ExperimentFolder = '/media/matthew/DataDisk/DrugKinetics_Data/20150813_DoxCompartmentModel_Take2'
ExperimentFolder = '/media/matthew/DataDisk/DrugKinetics_Data/20150902_DoxCompartmentModel_Take3';% good data for 231 uptake
%ExperimentFolder = '/media/matthew/DataDisk/DrugKinetics_Data/20151119_L12FluorescenceTest';

% Each experiment may have multiple data folders, put them all in
% DataFolder
%DataFolders{1} = fullfile(ExperimentFolder, ...
%    'Pathway_07292015','2015-07-29_002');
%DataFolders{2} = fullfile(ExperimentFolder, ...
%    'Pathway_07292015','2015-07-30_000');
%DataFolders{1} =  fullfile(ExperimentFolder, '2015-08-15_000');
DataFolders{1} =  fullfile(ExperimentFolder, '2015-09-02_000');
%DataFolders{1} =  fullfile(ExperimentFolder, '2015-11-19_000');

%PlateMapFile = fullfile(ExperimentFolder, ...
%    '20150727_MDAMB231_DoxorubicinCompartmentModel_PlateMap.csv');
%PlateMapFile = fullfile(ExperimentFolder, ...
%    '20150813_MDAMB231_DoxorubicinCompartmentModel_PlateMap.csv');
PlateMapFile = fullfile(ExperimentFolder, ...
    '20150902_MDAMB231_DoxorubicinCompartmentModel_PlateMap.csv');
%PlateMapFile = fullfile(ExperimentFolder, ...
%    '20151119_L12FluorescenceTest_PlateMap.csv');

% Read in plate map file to get concentrations and wells of interest
isDynamic = true;
[plateMap,DrugTime,cellType,folderList] = importPlateMap(PlateMapFile,isDynamic);

% Wells of interest...
% Define control wells
controlWells = cell2mat(cellfun(@(x) strcmp(x,'none'),cellType,'UniformOutput',false));
controlFolders = cell(sum(controlWells),1);
iter = 0;
for cwells = 1:length(controlWells)
    if controlWells(cwells)
        iter = iter+1;
        controlFolders{iter} = ['Well ' folderList{cwells}];
    end
end
%controlFolders = controlFolders(9:end);
controlDrugConcs = plateMap{3}(controlWells);

% Wells with cells
pbsWells = cell2mat(cellfun(@(x) strcmp(x,'null'),cellType,'UniformOutput',false));
cellWells = ~(pbsWells | controlWells);
cellFolders = cell(sum(cellWells),1);
iter = 0;
for cwells = 1:length(cellWells)
    if cellWells(cwells)
        iter = iter+1;
        cellFolders{iter} = ['Well ' folderList{cwells}];
    end
end
%cellFolders = cellFolders(1:8);
cellDrugConcs = plateMap{3}(cellWells);

% All folders of interest
foi = [controlFolders; cellFolders];
allConcs = cellDrugConcs;


%% Get time vectors and number of images

c_foi = foi{10};
[allTimes_hrs, numImage_perFolder] = createTimeVector(DataFolders,c_foi);

%% Read and process images

% fluorIntenMatrix = number of wells x number of images x 4
% 4: mean background, mean cell, mean control, std background
% need to include control and cell wells when making fluorIntenMatrix
fluorIntenMatrix = zeros(8,sum(numImage_perFolder),4);

% processing
% subtract off corresponding control well from cell well
% Signal is summation of 2 compartments (cytoplasm, nuclear)

maxInten = 2000;
fewImages = 19;
% Rows C and E
startDrug = 5-1;
endDrug = 28-1;
davec = 0:.01:1;

for controlWellIter = 1:10%sum(controlWells)
    disp(controlWellIter)
    
    for dfIter = 1:length(DataFolders)
        
        controlFolder_AIF = fullfile(DataFolders{dfIter},controlFolders{controlWellIter});
        controlFolder = fullfile(DataFolders{dfIter},controlFolders{controlWellIter+8});
        cellFolder = fullfile(DataFolders{dfIter},cellFolders{controlWellIter});
        %cellFolder
        %if ~exist(controlFolder,'dir')
        %    continue
        %end
        
        allOuttie =  zeros(numImage_perFolder(dfIter),4);
        
        for allIm = 0:fewImages%numImage_perFolder(dfIter)-1
            
            
            % correct images before processing
            flName = sprintf('Doxorubicin - n%06d.tif',allIm);
            aif_flName = fullfile(controlFolder_AIF,flName);
            %cell_flName = fullfile(cellFolder,flName);
            
            %Im_cell_current = double(imread(cell_flName))./maxInten;
            Im_aif = double(imread(aif_flName));
            
            %validCount = 0;
            %Im_cell = zeros(size(Im_cell_current));
            %Im_aif = zeros(size(Im_aif_current));
%             for maf = allIm-1:allIm+1
%                 
%                 % aif and cell wells correction
%                 flName = sprintf('Doxorubicin - n%06d.tif',maf);
%                 aif_flName = fullfile(controlFolder_AIF,flName);
%                 cell_flName = fullfile(cellFolder,flName);
%                 
%                 if (maf>=startDrug && maf<=endDrug) && ...
%                         (allIm>=startDrug && allIm<=endDrug)
%                     
%                     validCount = validCount+1;
%                     
%                     Im_cell_maf = double(imread(cell_flName))./maxInten;
%                     tpHist = hist(Im_cell_maf(:),davec);
%                     Im_cell = Im_cell + histeq(Im_cell_current,tpHist);
%                     
%                     
%                     Im_aif_maf = double(imread(aif_flName))./maxInten;
%                     tpHist = hist(Im_aif_maf(:),davec);
%                     Im_aif = Im_aif + histeq(Im_aif_current,tpHist);
%                     
%                     
%                 elseif (maf<startDrug || maf>endDrug) &&...
%                         (allIm<startDrug || allIm>endDrug) &&...
%                         (maf >=0)
%                     
%                     
%                     validCount = validCount+1;
%                     
%                     Im_cell_maf = double(imread(cell_flName))./maxInten;
%                     tpHist = hist(Im_cell_maf(:),davec);
%                     Im_cell = Im_cell + histeq(Im_cell_current,tpHist);
%                     
%                     
%                     Im_aif_maf = double(imread(aif_flName))./maxInten;
%                     tpHist = hist(Im_aif_maf(:),davec);
%                     Im_aif = Im_aif + histeq(Im_aif_current,tpHist);
%                     
%                     
%                 end
%             end
%             
%             Im_cell = Im_cell.*maxInten./validCount;
%             Im_aif = Im_aif.*maxInten./validCount;
            
            flNametl = sprintf('Transmitted Light - n%06d.tif',maf);
            cellTL_flName = fullfile(cellFolder,flNametl);
            Im_cellTL = double(imread(cellTL_flName));
            
%             flName = sprintf('Doxorubicin - n%06d.tif',allIm);
%             control_flName = fullfile(controlFolder,flName);
%             Im_control_current = double(imread(control_flName))./maxInten;
%             Im_control = zeros(size(Im_control_current));
%             validCount = 0;
%             if allIm < startDrug
%                 Im_control = Im_control_current.*maxInten;
%             else
%                 for maf = allIm-1:allIm+1
%                     if (maf>=startDrug) && (allIm>=startDrug)
%                         validCount = validCount + 1;
%                         flName = sprintf('Doxorubicin - n%06d.tif',maf);
%                         
%                         Im_control_maf = double(imread(cell_flName))./maxInten;
%                         tpHist = hist(Im_control_maf(:),davec);
%                         Im_control = Im_control + histeq(Im_control_current,tpHist);
%                     end
%                 end
%                 Im_control = Im_control.*maxInten./validCount;
%             end
            
            % take 3 images and process them!
            [outtie,SegOut] = processDynamicFluorescentImages(Im_aif,...
                Im_aif,Im_aif);
            %outtie(2) = mean(Im_aif(:));
            allOuttie(allIm+1,:) = outtie;
            
%             figure(13589);clf;
%             %subplot(131);imagesc(Im_cell);caxis([0 1000])
%             %imagesc(Im_cell-Im_aif);caxis([0 250])
%             %subplot(211)
%             %imagesc(Im_cellTL);caxis([1000 2000])
%             %colorbar
%             %daTitle = sprintf('%4.0f minutes',round((allTimes_hrs(allIm+1)-allTimes_hrs(5))*60))
%             %title(daTitle,'FontSize',16);
%             axis equal
%             axis off
%             
%             %subplot(212)
%             imagesc(Im_cell-Im_aif);caxis([0 250])
%             %subplot(132);imagesc(Im_control);caxis([0 1000]);title(num2str(allIm))
%             %subplot(133);imagesc(Im_aif);caxis([0 1000])
%             %colorbar
%             %daTitle = sprintf('%4.0f minutes',round((allTimes_hrs(allIm+1)-allTimes_hrs(5))*60))
%             %title(daTitle,'FontSize',16);
%             set(gcf,'color','w');
%             axis equal
%             axis off
%             
%             frame = getframe(13589);
%             im = frame2im(frame);
%             [imind,damap] = rgb2ind(im,256);
% 
%         if allIm == 3;
%               imwrite(imind,damap,'DoxUptake.gif','gif', 'Loopcount',inf, 'DelayTime',.75);
%         else
%               imwrite(imind,damap,'DoxUptake.gif','gif','WriteMode','append', 'DelayTime',.75);
%         end
%             
            %pause(.01)
            
        end
        
        startInsert = sum(numImage_perFolder(1:dfIter-1)) + 1;
        endInsert = sum(numImage_perFolder(1:dfIter));
        fluorIntenMatrix(controlWellIter, startInsert:endInsert, :) = allOuttie;
        
    end
    
end

%%
% look at fluorescent signal over time
clrs = {'r','g','b','k','r','g','b','k','r','g','b','k'};
figure(997);clf
for aaa = 1:8
    %figure(997);clf
    hold on
    %plot(allTimes_hrs, fluorIntenMatrix(aaa,:,1),'g.')
    %plot(allTimes_hrs, fluorIntenMatrix(aaa,:,2),'ro')
    plot(allTimes_hrs, fluorIntenMatrix(aaa,:,3),clrs{aaa})
    xlabel('Time (hrs)')
    ylabel('Intensity')
    axis([0 10 0 2000])
    %legend('Intracellular','Extracellular','Control')
    hold off
    
    %pause(1)
end

legend(num2str(controlDrugConcs(1:8)))

%%

% convert to concentration prior to fitting model parameters
fluorIntenMatrix_mod = fluorIntenMatrix(:,1:fewImages+1,:);
%fluorIntenMatrix_mod(:, :, 2) = fluorIntenMatrix_mod(:,:, 3);
concMatrix = convertToConcentration(fluorIntenMatrix_mod, allConcs, allTimes_hrs);

concMatrix(:,:,1) = concMatrix(:,:,1) - concMatrix(:,:,2);
concMatrix(:,29:end,2) = 0;

allTimes_mod = allTimes_hrs(1:fewImages+1);
%concMatrix(:,1:4,:) = 0;
%concMatrix(:,25:28,:) = [];
%allTimes_mod(25:28) = [];

%concMatrix(:,:,1) = concMatrix(:,:,1)-concMatrix(:,:,2);

figure(997);clf
figure(998);clf
for aaa = 1:4%1:4
    controlDrugConcs(aaa)
    options=optimset('TolFun',1e-6,'TolX',1e-6,'MaxIter',1e50,'Display','off','PrecondBandWidth',0);
    
    % Add in weighting vector?
    weight_v = ones(1,size(concMatrix,2));
    
    %fluorIntenMatrix(aaa,1:28,1) = 0;
    %fluorIntenMatrix(aaa,1:28,2) = 0;
    
    
    % 2 compartment 3 parameter
    % kpars =   [k1 k2 k3]
    eod = 0;
    estim =     [2 0.05 .05];
    lowBound =  [0.001  0.001 .001];
    highBound = [10  5 5];
    [parm,Rn,Res,whyExit,~,~,J] = lsqcurvefit(@(kpars,cp) modelCt3C3P(kpars,cp,allTimes_mod(1:end-eod)),...
        estim, squeeze(concMatrix(aaa,1:end-eod,2)), squeeze(concMatrix(aaa,1:end-eod,1)),...
        lowBound, highBound, options);
    ks = parm;
    ks(2) = parm(1)/parm(2);
    
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
    
    figure(998)
    subplot(131)
    semilogx(0,0)
    hold on
    errorbar(allConcs(aaa),z_ci(1,2),tst*se_g(1),'ro','MarkerSize',10,'MarkerFaceColor',[1 0 0],'LineWidth',2)
    
    subplot(132)
    semilogx(0,0)
    hold on
    errorbar(allConcs(aaa),z_ci(2,2),tst*se_g(2),'ko','MarkerSize',10,'MarkerFaceColor',[0 0 0],'LineWidth',2)
    
    subplot(133)
    semilogx(0,0)
    hold on
    errorbar(allConcs(aaa),z_ci(3,2),tst*se_g(3),'bo','MarkerSize',10,'MarkerFaceColor',[0 0 1],'LineWidth',2)
    
    %[parm' tst*se_g]
    %se_g
    
    %disp([parm allConcs(aaa)])modelCt2C3P
    
    fwdModEval = modelCt3C3P(parm,squeeze(concMatrix(aaa,1:end-eod,2)),allTimes_mod(1:end-eod));
    figure(997)
    if aaa>0 && aaa<2
    %subplot(2,2,aaa)
    hold on
    plot(allTimes_mod(1:end-eod), squeeze(concMatrix(aaa,1:end-eod,1)),'ko','MarkerSize',10,'MarkerFaceColor',[0 0 0],'LineWidth',2)
    plot(allTimes_mod(1:end-eod), squeeze(concMatrix(aaa,1:end-eod,2)),'bo','MarkerSize',10,'MarkerFaceColor',[0 0 1],'LineWidth',2)
    %plot(allTimes_hrs(1:end-eod), squeeze(fluorIntenMatrix(aaa,1:end-eod,3))-1000,'ks')
    plot(allTimes_mod(1:end-eod), fwdModEval, 'r-','LineWidth',5)
    %title([num2str(controlDrugConcs(aaa)) ' nM'])
    legend('Intracellular','Extracellular','Model Fit','Location','Best')
    xlabel('Time (hrs)','FontSize',32)
    ylabel('Concentration (nM)','FontSize',32)
    set(gca,'FontSize',22)
    axis([-.1 20 -100 3000*2^-(aaa-1)])
    hold off
    end
    
    %pause
    
end
figure(998)
subplot(131)
%xlabel('Concentration (nM)','FontSize',24)
ylabel('k_{12}','FontSize',24)
axis([0 2600 0 .1])
set(gca,'FontSize',15)

subplot(132)
xlabel('Concentration (nM)','FontSize',24)
ylabel('V_{d}','FontSize',24)
axis([0 2600 0 .2])
set(gca,'FontSize',15)

subplot(133)
%xlabel('Concentration (nM)','FontSize',24)
ylabel('k_{23}','FontSize',24)
axis([0 2600 0 .1])
set(gca,'FontSize',15)

% % % %%
% % % % plot change in dox signal over time
% % %
% % % figure(1);clf;
% % % hold on
% % % clrs = {'r.-','gs-','bo-',...
% % %     'rd-','g.-','bs-',...
% % %     'ro-','gd-'};
% % % for a = 1:size(fluorIntenMatrix,1)
% % %
% % %     %fluorIntenMatrix(a,1:28,2) = 0;
% % %
% % %     plot(allTimes_hrs,fluorIntenMatrix(a,:,2),clrs{mod(a,8)+1})
% % % end
% % % xlabel('Times (hr)')
% % % ylabel('Median Fluorescence')
% % % title('Doxorubicin Fluorescence, Background')
% % % hold off
% % %
% % % %%
% % %
% % % % fit linear model to fluorescence vs dox concentration data
% % % figure(2);clf;
% % % allConcs = controlDrugConcs;
% % %
% % % par = polyfit(allConcs(2:8),fluorIntenMatrix(2:8,30,2),1);
% % % cMod = linspace(0,10000);
% % % fMod = par(1)*cMod + par(2);
% % %
% % % hold on
% % % plot(allConcs(1:8), fluorIntenMatrix(1:8,30,2),'ro')
% % % plot(cMod,fMod,'b-')
% % % errorbar(allConcs(1:8),fluorIntenMatrix(1:8,30,2),fluorIntenMatrix(1:8,30,4),'r.')
% % % xlabel('Concentration (nM)')
% % % ylabel('Fluorescence Signal')
% % % axis([-200 10200 0 4500])
% % %
% % % %% Segment cells using TL, use mask to get timecourse for doxorubicin
% % %
% % % % outputs:
% % % % average background signal (extracellular space)
% % % % average uptake within all cells in FOV
% % %
% % % for allIm = 0:numImages-1
% % %
% % %     %flName = sprintf('./CorrectedIm/Doxorubicin - n0000%02d.tif',allIm);
% % %     %flName = sprintf('./Well C03/Doxorubicin - n%06d.tif',allIm);
% % %     flName = sprintf('./Well C03/Transmitted Light - n%06d.tif',allIm);
% % %
% % %     flName = sprintf('./Well C03/Transmitted Light - n%06d.tif',allIm);
% % %
% % %     Im_c = imread(flName);
% % %
% % %     [~, threshold] = edge(Im_c, 'sobel');
% % %     fudgeFactor = .75;
% % %     BWs = edge(Im_c,'sobel', threshold * fudgeFactor);
% % %
% % %     se90 = strel('line', 3, 90);
% % %     se0 = strel('line', 3, 0);
% % %     BWsdil = imdilate(BWs, [se90 se0]);
% % %     BWdfill = imfill(BWsdil, 'holes');
% % %
% % %     seD = strel('diamond',2);
% % %     BWdfill = imerode(BWdfill,seD);
% % %
% % %     BWnobord = imclearborder(BWdfill, 4);
% % %
% % %     seD = strel('diamond',1);
% % %     BWfinal = imerode(BWnobord,seD);
% % %     BWfinal = imerode(BWfinal,seD);
% % %
% % %     cc = bwconncomp(BWfinal);
% % %     numPixels = cellfun(@numel,cc.PixelIdxList);
% % %     noise = find(numPixels < 50);
% % %     for aa = 1:length(noise)
% % %         BWfinal(cc.PixelIdxList{noise(aa)}) = 0;
% % %     end
% % %
% % %
% % %     BWoutline = bwperim(BWfinal);
% % %     Segout = Im_c;
% % %     Segout(BWoutline) = max(Im_c(:));
% % %
% % %     figure(1);clf;
% % %     subplot(121)
% % %     imagesc(Im_c)
% % %     subplot(122)
% % %     imagesc(Segout)
% % %
% % %     pause(1)
% % %
% % % end
% % %
% % %
% % % %% Using timecourses, fit model to data
% % %
% % % % fitDoxCompartmentModel(backgroundSignal, cellSignal)
