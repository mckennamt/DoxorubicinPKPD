


% Dynamic dox modeling
close all;clear all;clc

CellLine = 'SUM149';
currDate = 20161022;
SelectCompartmentModelExperiment

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

c_foi = foi{1};
[allTimes_hrs, numImage_perFolder] = createTimeVector(DataFolders,c_foi);
allTimes_hrs = allTimes_hrs(1:end-1);
%numImage_perFolder = numImage_perFolder;

%% Read and process images, single exposure time

%preprocessedFluor = load(fullfile(pwd, ['SegmentationMasks_20160803_' CellLine '_headless.mat']),'CellMasks','allSeeds');
%CellMasks_fluorSeg = preprocessedFluor.CellMasks;
%allSeeds_fluorSeg = preprocessedFluor.allSeeds;

CellMasks = cell(length(cellFolders),sum(numImage_perFolder));
allSeeds = cell(length(cellFolders),sum(numImage_perFolder));

%load(fullfile(pwd, ['SegmentationMasks_' CellLine '_complete.mat']));
%load(fullfile(pwd, ['SegmentationMasks_' CellLine '_20160914_complete.mat']));
%load(fullfile(pwd, ['SegmentationMasks_' CellLine '_' num2str(currDate) '_complete.mat']));
numCells = 50;

%preprocessedFluor = load(fullfile(pwd, ['SegmentationMasks_20160803_' CellLine '_headless.mat']),'CellMasks','allSeeds');
%preprocessedFluor = load(fullfile(pwd, ['SegmentationMasks_20160914_' CellLine '_headless.mat']),'CellMasks','allSeeds');
preprocessedFluor = load(fullfile(pwd, ['SegmentationMasks_' num2str(currDate) '_' CellLine '_headless.mat']),'CellMasks','allSeeds');

CellMasks_fluorSeg = preprocessedFluor.CellMasks;
allSeeds_fluorSeg = preprocessedFluor.allSeeds;

for wellIter = 1:length(cellFolders)
    disp(wellIter)
    
    for dfIter = 1:length(DataFolders)
        
        cellFolder = fullfile(DataFolders{dfIter},cellFolders{wellIter});
        aifFolder = fullfile(DataFolders{dfIter},aifFolders{wellIter});
        previousSeeds = [];
        previousTL = [];
        croppedCell = [];
        
        %for allIm = 0:numImage_perFolder(dfIter)-1
        for allIm = 0:60
            disp(allIm)
            
            %if wellIter<7
            %    continue
            %elseif  wellIter == 7 && allIm==54%~isempty(CellMasks{wellIter,allIm+1})
            %    load(fullfile(pwd, ['SegmentationMasks_' CellLine '_complete.mat']),'previousSeeds');
            %    continue
            %elseif wellIter == 7 && allIm < 55
            %    continue
            %end
            if ~isempty(CellMasks{wellIter,allIm+1})
                continue
            end
            
            if allConcs(wellIter) == 0
                CellMasks{wellIter,allIm+1} = CellMasks_fluorSeg{wellIter,allIm+1};
                allSeeds{wellIter,allIm+1} = allSeeds_fluorSeg{wellIter,allIm+1};
                continue
            end
            
            
            
            flNametl = sprintf('Transmitted Light - n%06d.tif',allIm);
            cellTL_flName = fullfile(cellFolder,flNametl);
            Im_cellTL = double(imread(cellTL_flName));
            
            flNametl = sprintf('Doxorubicin - n%06d.tif',allIm);
            cellTL_flName = fullfile(cellFolder,flNametl);
            Im_cellFL = double(imread(cellTL_flName));
            
            flNametl = sprintf('Doxorubicin - n%06d.tif',allIm);
            cellTL_flName = fullfile(aifFolder,flNametl);
            Im_aifFL = double(imread(cellTL_flName));
            
            probGoodFluorescentSignal = allIm >= (startDrug) & cellDrugConcs(wellIter)>0;
            
            [cellMask,currentSeeds, croppedCell] = segmentCells_main(Im_cellTL,...
                Im_cellFL, Im_aifFL, probGoodFluorescentSignal, ...
                numCells, previousSeeds, previousTL,...
                CellMasks_fluorSeg{wellIter,allIm+1},...
                allSeeds_fluorSeg{wellIter,allIm+1},croppedCell);
            
            numCells = size(currentSeeds,1);
            previousSeeds = currentSeeds;
            previousTL = Im_cellTL;
            
            CellMasks{wellIter,allIm+1} = cellMask;
            allSeeds{wellIter,allIm+1} = currentSeeds;
            
        end
        
    end
    %save(fullfile(pwd, ['SegmentationMasks_' CellLine '_complete.mat']));
    
end

%save(fullfile(pwd, ['SegmentationMasks_' CellLine '_20160914_complete.mat']),'-v7.3');
save(fullfile(pwd, ['SegmentationMasks_' CellLine '_' num2str(currDate) '_complete.mat']), '-v7.3');

backupDir = '/home/matthew/Documents/Vanderbilt/Gs/model code and data';
save(fullfile(backupDir, ['SegmentationMasks_' CellLine '_' num2str(currDate) '_complete.mat']), '-v7.3');


% %% run quick check of uptake signal, just using masks
%
% daSignal = [];
%
% for wellIter = 1%:length(cellFolders)
%     disp(wellIter)
%
%     for dfIter = 1:length(DataFolders)
%
%         cellFolder = fullfile(DataFolders{dfIter},cellFolders{wellIter});
%         aifFolder = fullfile(DataFolders{dfIter},aifFolders{wellIter});
%
%         for allIm = 0:numImage_perFolder(dfIter)-1
%
%             flNametl = sprintf('Doxorubicin - n%06d.tif',allIm);
%             cellTL_flName = fullfile(cellFolder,flNametl);
%             Im_cellFL = double(imread(cellTL_flName));
%
%             flNametl = sprintf('Doxorubicin - n%06d.tif',allIm);
%             cellTL_flName = fullfile(aifFolder,flNametl);
%             Im_aifFL = double(imread(cellTL_flName));
%
%             diffIm = Im_cellFL - Im_aifFL;
%
%             daSignal(allIm+1,1) = mean(diffIm(CellMasks{wellIter,allIm+1}));
%             daSignal(allIm+1,2) = mean(Im_aifFL(CellMasks{wellIter,allIm+1}));
%
%             %figure(1);clf;
%             %imagesc(diffIm)
%             %caxis([-20 100])
%             %hold on
%             %h = imshow(cat(3,ones(size(Im_cellTL)),zeros(size(Im_cellTL)),zeros(size(Im_cellTL))));
%             %set(h,'AlphaData',.2.*CellMasks{wellIter,allIm+1});
%             %pause
%
%         end
%
%     end
%
% end
% %%
% figure(2);clf;
% hold on
% plot(allTimes_hrs(:)'-allTimes_hrs(startDrug),daSignal(:,1),'g.-')
% plot(allTimes_hrs(:)'-allTimes_hrs(startDrug),daSignal(:,2)-daSignal(1,2),'r.-')
% legend('Cell Uptake','Dox AIF','Location','Best')
