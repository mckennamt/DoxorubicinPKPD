


% Dynamic dox modeling
close all;clear all;clc

CellLine = 'MDAMB468';
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
numImage_perFolder = numImage_perFolder-1;

%% Read and process images, single exposure time

%preprocessedFluor = load(fullfile(pwd, ['SegmentationMasks_20160803_' CellLine '_headless.mat']),'CellMasks','allSeeds');
%CellMasks_fluorSeg = preprocessedFluor.CellMasks;
%allSeeds_fluorSeg = preprocessedFluor.allSeeds;

CellMasks = cell(length(cellFolders),sum(numImage_perFolder));
allSeeds = cell(length(cellFolders),sum(numImage_perFolder));

load(fullfile(pwd, ['SegmentationMasks_' CellLine '_complete.mat']));
%load(fullfile(pwd, ['SegmentationMasks_' CellLine '_20160914_complete.mat']));
numCells = 40;

%preprocessedFluor = load(fullfile(pwd, ['SegmentationMasks_20160803_' CellLine '_headless.mat']),'CellMasks','allSeeds');
%preprocessedFluor = load(fullfile(pwd, ['SegmentationMasks_20160914_' CellLine '_headless.mat']),'CellMasks','allSeeds');
CellMasks_fluorSeg = preprocessedFluor.CellMasks;
allSeeds_fluorSeg = preprocessedFluor.allSeeds;


%%
for wellIter = 1
    disp(wellIter)
    
    for dfIter = 1
        
        cellFolder = fullfile(DataFolders{dfIter},cellFolders{wellIter});
        aifFolder = fullfile(DataFolders{dfIter},aifFolders{wellIter});
        previousSeeds = [];
        previousTL = [];
        croppedCell = [];
        
        for allIm = 30%0:numImage_perFolder(dfIter)-1
            disp(allIm)
            
            
            
            flNametl = sprintf('Transmitted Light - n%06d.tif',allIm);
            cellTL_flName = fullfile(cellFolder,flNametl);
            Im_cellTL = double(imread(cellTL_flName));
            
            flNametl = sprintf('Doxorubicin - n%06d.tif',allIm);
            cellTL_flName = fullfile(cellFolder,flNametl);
            Im_cellFL = double(imread(cellTL_flName));
            
            flNametl = sprintf('Doxorubicin - n%06d.tif',allIm);
            cellTL_flName = fullfile(aifFolder,flNametl);
            Im_aifFL = double(imread(cellTL_flName));
            
            alVal = .2;
            
            f2 = figure(2);clf;
            imagesc(Im_cellTL)
            hold on
            lab = label2rgb(CellMasks{wellIter,allIm+1},'autumn');
            caxis([1000 1500])
            h = imshow(lab);
            set(h,'AlphaData',alVal)
            colormap gray
            set(f2,'Position',[80 1080 570 425])
            
            f3 = figure(3);clf;
            imagesc(Im_cellFL)
            colormap gray
            %caxis([240 400])
            axis off
            hold on
            h = imshow(lab);
            set(h,'AlphaData',alVal)
            set(f3,'Position',[95 565 570 425])
            
            
            
        end
        
    end
    
end



