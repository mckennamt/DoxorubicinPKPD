


% process all fluorescent images, no user input

% Dynamic dox modeling
close all;clear all;clc

%allCellLines = {'SUM149','MDAMB453','MDAMB468','MDAMB231'};
allCellLines = {'SUM149'};
%allDates = [20161018, 20161022];
allDates = [20161022];

for cellLineIter = 1:length(allCellLines)
    for dateIter = 1:length(allDates)
        currDate = allDates(dateIter);
        
        CellLine = allCellLines{cellLineIter};
        
        clearvars -except allCellLines CellLine cellLineIter currDate allDates dateIter
        
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
        allTimes_hrs = allTimes_hrs - allTimes_hrs(startDrug+1);
        %numImage_perFolder = numImage_perFolder-1;
        
        %% Read and process images, single exposure time
        
        CellMasks = cell(length(cellFolders),sum(numImage_perFolder));
        allSeeds = cell(length(cellFolders),sum(numImage_perFolder));
        numCells = 40;
        fprintf('%s\n',CellLine)
        for wellIter = 1:length(cellFolders)
            
            fprintf('Well iter:\t%2d\n',wellIter)
            
            for dfIter = 1:length(DataFolders)
                
                cellFolder = fullfile(DataFolders{dfIter},cellFolders{wellIter});
                aifFolder = fullfile(DataFolders{dfIter},aifFolders{wellIter});
                
                parfor allIm = 0:numImage_perFolder(dfIter)-1
                %for allIm = 0:numImage_perFolder(dfIter)-1
                    
                    %fprintf('%2d ',allIm)
                    
                    flNametl = sprintf('Doxorubicin - n%06d.tif',allIm);
                    cellTL_flName = fullfile(cellFolder,flNametl);
                    Im_cellFL = double(imread(cellTL_flName));
                    
                    flNametl = sprintf('Doxorubicin - n%06d.tif',allIm);
                    cellTL_flName = fullfile(aifFolder,flNametl);
                    Im_aifFL = double(imread(cellTL_flName));
                    
                    if allConcs(wellIter) == 0
                        cellMask = false(size(Im_cellFL));
                        currentSeeds = [500 500];
                    else
                        [cellMask,currentSeeds] = segmentFluorescentDoxImages_headless(Im_cellFL,Im_aifFL);
                    end
                    
                    CellMasks{wellIter,allIm+1} = cellMask;
                    allSeeds{wellIter,allIm+1} = currentSeeds;
                    
                end
                
            end
            fprintf('\n')
            
        end
        
        fprintf('\n\n')
        
        save(fullfile(pwd, ['SegmentationMasks_' num2str(currDate) '_' CellLine '_headless.mat']), '-v7.3');
        
    end
end