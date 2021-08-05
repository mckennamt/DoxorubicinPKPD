
% Dynamic dox modeling
close all;clear all;clc

%CellLine = 'SUM149';
%currDate = 20161022;
%SelectCompartmentModelExperiment

%load(fullfile(pwd, ['SegmentationMasks_' CellLine '_complete.mat']));
%load(fullfile(pwd, ['SegmentationMasks_' CellLine '_20160914_complete.mat']));
%load(fullfile(pwd, ['SegmentationMasks_' CellLine '_' num2str(currDate) '_complete.mat']));

%SelectCompartmentModelExperiment

allCellLines = {'SKMEL5'};
allDates = [20161110];
ExperimentFolder = '/media/matthew/DataDisk/DrugKinetics_Data/20161110_SKMEL5_SC07_Tariquidar';
DataFolders{1} =  fullfile(ExperimentFolder, '2016-11-10_000');
PlateMapFile = fullfile(ExperimentFolder, ...
    '20161110_SKMEL5_SC07_Tariquidar_PlateMap.csv');
startDrug = 3;%Image file number (first image with drug on cells; zero-index)
endDrug = [25 25];%Image file number (last image with drug on cells; zero-index)

currDate = allDates(1);
CellLine = allCellLines{1};

load(fullfile(pwd, ['SegmentationMasks_' num2str(currDate) '_' CellLine '_H2BRFP.mat']));

%%

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

analysisWells = [1:9 11:19];
%analysisWells = [1:8 11:18];

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

%maxInten = 1;%1000
%davec = 0:.01:1; %used for histogram correction


for wellIter = analysisWells
    fprintf('%2d ',wellIter)
    
    for dfIter = 1:length(DataFolders)
        
        controlFolder_blank = fullfile(DataFolders{dfIter},aifFolders{9});
        controlFolder_blank2 = fullfile(DataFolders{dfIter},aifFolders{10});
        
        controlFolder_AIF = fullfile(DataFolders{dfIter},aifFolders{wellIter});
        cellFolder = fullfile(DataFolders{dfIter},cellFolders{wellIter});
        drugOnlyFolder = fullfile(DataFolders{dfIter},drugOnlyFolders{wellIter});
        
        allOuttie =  zeros(numImage_perFolder(dfIter),4);
        fewImages = numImage_perFolder(dfIter)-1;
        
        for allIm = 0:fewImages
            
            % Need to re-write processDynamicFluorescentImages to do the
            % histogram normalization technique
            % correct images before processing
            
            flName = sprintf('Doxorubicin - n%06d.tif',allIm);
            
            aif_blank_flName = fullfile(controlFolder_blank,flName);
            aif_blank2_flName = fullfile(controlFolder_blank2,flName);
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
            %             %subplot(121)
            %imagesc(Im_cell);caxis([100 800]);
            %             %subplot(122)
            %             %imagesc(Im_cell - Im_aif);caxis([-100 100]);
            %colormap gray
            %hold on
            %lab = label2rgb(cCellMask,'autumn');
            %h = imshow(lab);
            %set(h,'AlphaData',.1)
            %             %subplot(122)
            %             %imDiff = Im_cell - Im_aif;
            %             %%imDiff(imDiff<-20) = 0;
            %             %%imDiff(imDiff>100) = 0;
            %             %[n,x] = hist(imDiff(cCellMask),-20:10:100);
            %             %plot(clip_mean, linspace(0,1,1000),'r-','LineWidth',5)
            %             %hold on
            %             %bar(x,n./max(n),'g')
            %             %axis([-20 100 0 1.1])
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
fprintf('\n\n')

%% for each timepoint, make concentration look up table

concFluorInten_Map = zeros(length(allTimes_hrs),2);
concMatrix = zeros(size(fluorIntenMatrix));

for dfIter = 1:length(DataFolders)
    fewImages = numImage_perFolder(dfIter)-1;
    
    for allIm = 0:fewImages
        fprintf('%2d ', allIm)
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
        
        % convert intensity to concentration
        if allIm > startDrug
            concMatrix(:,allIm+1,:) = ...
                (fluorIntenMatrix(:,allIm+1,:) - linFit(2))./linFit(1);
        elseif allIm == startDrug
            for cleanUp = 0:(startDrug)
                concMatrix(:,cleanUp+1,:) = ...
                    (fluorIntenMatrix(:,cleanUp+1,:) - linFit(2))./linFit(1);
            end
        end
        
        %lb = -1.96*stdInten;
        %ub = 1.96*stdInten;
        %cf = figure(1);clf;
        %figure(121);clf
        %plot(uConc,fluorInten,'r.')
        %pause
        %ss = [1 4:length(uConc)];
        %hold on
        %cMod = linspace(0,2800);
        %fMod = linFit(1)*cMod + linFit(2);
        %plot(cMod,fMod,'k-','LineWidth',1.5);
        %errorbar(uConc(ss), fluorInten(ss), lb(ss), ub(ss),'ro')
        %axis([-100 3000 200 800]);
        %%set(cf,'Position',[100 100 900 250])
        %set(gca,'FontSize',18)
        %pause(.5)
    end
    
end

fprintf('\n\n')
%%
if strcmp(CellLine,'MDAMB453')
    concMatrix(:,54,:) = [];
    fluorIntenMatrix(:,54,:) = [];
    allTimes_hrs(54) = [];
end

if strcmp(CellLine,'MDAMB468')
    concMatrix(:,end,:) = [];
    fluorIntenMatrix(:,end,:) = [];
    allTimes_hrs(end) = [];
end

save(['IntensityTimecourses_' num2str(currDate) '_' CellLine '.mat'],'concMatrix',...
    'allTimes_hrs','analysisWells',...
    'startDrug','endDrug','aifDuration','drugOnlyConcs',...
    'allConcs','concFluorInten_Map','fluorIntenMatrix')


