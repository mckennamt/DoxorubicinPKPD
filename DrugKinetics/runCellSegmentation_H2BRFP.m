


% process all fluorescent images, no user input

% Dynamic dox modeling
close all;clear all;clc

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

%SelectCompartmentModelExperiment

%% Read in plate map file to get concentrations and wells of interest
isDynamic = true;
%[plateMap,DrugTime,cellType,folderList] = importPlateMap(PlateMapFile,isDynamic);
[plateMap,DrugTime,cellType,folderList] = importPlateMap_Sensitizer(PlateMapFile,isDynamic,true);

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
%allTimes_hrs = allTimes_hrs(1:end-1);
%numImage_perFolder = numImage_perFolder-1;

%% Read and process images, single exposure time

CellMasks = cell(length(cellFolders),sum(numImage_perFolder));
allSeeds = cell(length(cellFolders),sum(numImage_perFolder));

fprintf('%s\n',CellLine)
for wellIter = 1:length(cellFolders)
    
    fprintf('Well iter:\t%2d\n',wellIter)
    
    for dfIter = 1:length(DataFolders)
        
        cellFolder = fullfile(DataFolders{dfIter},cellFolders{wellIter});
        aifFolder = fullfile(DataFolders{dfIter},aifFolders{wellIter});
        
        parfor allIm = 0:numImage_perFolder(dfIter)-1
            %for allIm = 0:numImage_perFolder(dfIter)-1
            
            flNametl = sprintf('DsRed - n%06d.tif',allIm);
            cellTL_flName = fullfile(cellFolder,flNametl);
            Im_cellH2B = double(imread(cellTL_flName));
            [cellMask,currentSeeds] = segment_H2BRFPImages(Im_cellH2B);
            
            CellMasks{wellIter,allIm+1} = cellMask;
            allSeeds{wellIter,allIm+1} = currentSeeds;
        end
        
    end
    fprintf('\n')
    
end

fprintf('\n\n')

save(fullfile(pwd, ['SegmentationMasks_' num2str(currDate) '_' CellLine '_H2BRFP.mat']), '-v7.3');
