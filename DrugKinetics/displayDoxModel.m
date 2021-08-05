
% Dynamic dox modeling
close all;clear all;clc

CellLine = 'MDAMB231';
SelectCompartmentModelExperiment

load(fullfile(pwd, ['SegmentationMasks_' CellLine '_complete.mat']));

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


%% Get time vectors and number of images

% Should be loaded with cell segmentation data
c_foi = foi{1};
[allTimes_hrs, numImage_perFolder] = createTimeVector(DataFolders,c_foi);
allTimes_hrs = allTimes_hrs - allTimes_hrs(startDrug+1);

%% Read and process images, single exposure time

% fluorIntenMatrix = number of wells x number of images x 4
% 4: mean cell value, mean in AIF, std cell, std AIF
fluorIntenMatrix = zeros(1,sum(numImage_perFolder),4);
%fluorIntenMatrix = zeros(8,50,4);

% processing
% subtract off corresponding aif well from cell well
% Signal is summation of 2 compartments (cytoplasm, nuclear)

maxInten = 1;%1000
davec = 0:.01:1; %used for histogram correction


for wellIter = 1%1:20
    disp(wellIter)
    
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
            
            allOuttie(allIm+1,:) = [mean(Im_cell(cCellMask)) mean(Im_cell(~cCellMask)),...
                mean(Im_aif(:)) mean(Im_blank(:))];
            
            %cf = figure(13589);clf;
            %imagesc(Im_cell - Im_aif);caxis([-10 100])
            %colormap gray
            %set(cf,'Position',[1200 200 500 400])
            %set(gca,'visible','off');
            %disp(allTimes_hrs(allIm+1))
            %pause
            
            %OLD STUFF, used to make gifs
            %cf = figure(13589);clf;
            %subplot(121)
            %daTitle = sprintf('%4.2f hours',round((allTimes_hrs(allIm+1)-allTimes_hrs(startDrug+1)),2));
            %imagesc(Im_cell - Im_aif);caxis([-10 100]);title(daTitle,'FontSize',8);
            %colormap gray
            %hold on
            %lab = label2rgb(cCellMask,'autumn');
            %h = imshow(lab);
            %set(h,'AlphaData',.1)
            %subplot(122)
            %[n,x] = hist(Im_cell(cCellMask)-mean(Im_aif(:)),-10:10:100);
            %bar(x,n./sum(n))
            %axis([-20 100 0 1])
            %%subplot(131);imagesc(Im_cell-Im_aif);caxis([-20 200]);title(daTitle,'FontSize',8);
            %%subplot(132);imagesc(Im_aif);caxis([200 600]);title(daTitle,'FontSize',8);
            %%subplot(133);imagesc(Im_cellTL);caxis([1000 2000]);title(daTitle,'FontSize',8);
            %
            %subplot(2,1,1);imagesc(Im_cellTL);caxis([1000 2000]);title(daTitle,'FontSize',8);
            %axis off
            %subplot(2,1,2);imagesc(Im_cell-Im_aif);caxis([-20 100]);
            %axis off
            %colormap gray
            %%colorbar
            %%daTitle = sprintf('%4.0f minutes',round((allTimes_hrs(allIm+1)-allTimes_hrs(5))*60));
            %%title(daTitle,'FontSize',16);
            %%set(gcf,'color','w');
            %%axis equal
            %%axis off
            %
            %%subplot(212)
            %imagesc(Im_cell-Im_aif);caxis([0 250])
            %%subplot(132);imagesc(Im_control);caxis([0 1000]);title(num2str(allIm))
            %%subplot(133);imagesc(Im_aif);caxis([0 1000])
            %%colorbar
            %%daTitle = sprintf('%4.0f minutes',round((allTimes_hrs(allIm+1)-allTimes_hrs(5))*60))
            %%title(daTitle,'FontSize',16);
            %
            %%set(cf,'Position',[100 100 1800 600])
            %set(cf,'Position',[100 100 400 600])
            %set(gca,'visible','off');
            %
            %frame = getframe(cf);
            %im = frame2im(frame);
            %[imind,damap] = rgb2ind(im,256);
            %
            %if allIm == 0;
            %    imwrite(imind,damap,'DoxUptake_MDAMB468_2500_CQSB.gif','gif', 'Loopcount',inf, 'DelayTime',.75);
            %else
            %    imwrite(imind,damap,'DoxUptake_MDAMB468_2500_CQSB.gif','gif','WriteMode','append', 'DelayTime',.75);
            %end
            %
            %pause%(.001)
            
            
        end
        
        startInsert = sum(numImage_perFolder(1:dfIter-1)) + 1;
        endInsert = sum(numImage_perFolder(1:dfIter));
        fluorIntenMatrix(wellIter, startInsert:endInsert, :) = allOuttie;
        
    end
    
    tc_diff = fluorIntenMatrix(wellIter,:,1)-fluorIntenMatrix(wellIter,:,2);
    min_tcdiff = abs(min(tc_diff));
    
    figure(997);clf
    hold on
    plot(allTimes_hrs, fluorIntenMatrix(wellIter,:,1),'g.')
    plot(allTimes_hrs, fluorIntenMatrix(wellIter,:,2),'b.')
    plot(allTimes_hrs, tc_diff,'ro')
    xlabel('Time (hrs)')
    ylabel('Intensity')
    %axis([0 24 200 500])
    %axis([0 24 -10 30])
    %legend('Intracellular','Extracellular','Control')
    hold off
    
end

concMatrix = fluorIntenMatrix;
%%
cf = figure(1212);clf;
plot(allTimes_hrs,tc_diff,'ro')
axis([-2.6 22 -3 40])
%xlabel('Time (hrs)')
%ylabel('Intensity')
set(gca,'FontSize',12)
set(cf,'Position',[1200 200 1000 200])



