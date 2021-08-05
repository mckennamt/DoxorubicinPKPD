
% Dynamic dox modeling
%close all;clear all;clc

% Define experiment folder
%ExperimentFolder = '/media/matthew/DataDisk/DrugKinetics_Data/20150727_DoxCompartmentModel_Take1';
%ExperimentFolder = '/media/matthew/DataDisk/DrugKinetics_Data/20150813_DoxCompartmentModel_Take2'
%ExperimentFolder = '/media/matthew/DataDisk/DrugKinetics_Data/20150902_DoxCompartmentModel_Take3';% good data for 231 uptake
%ExperimentFolder = '/media/matthew/DataDisk/DrugKinetics_Data/20151119_L12FluorescenceTest';
ExperimentFolder = '/media/matthew/DataDisk/DrugKinetics_Data/20160315_FixDoxSignalExperiment';


% Each experiment may have multiple data folders, put them all in
% DataFolder
%DataFolders{1} = fullfile(ExperimentFolder, ...
%    'Pathway_07292015','2015-07-29_002');
%DataFolders{2} = fullfile(ExperimentFolder, ...
%    'Pathway_07292015','2015-07-30_000');
%DataFolders{1} =  fullfile(ExperimentFolder, '2015-08-15_000');
%DataFolders{1} =  fullfile(ExperimentFolder, '2015-09-02_000');

%DataFolders{1} =  fullfile(ExperimentFolder, '2015-11-19_000');
%DataFolders{1} =  fullfile(ExperimentFolder, '2015-11-19_001');
DataFolders{1} =  fullfile(ExperimentFolder, '2016-03-15_003');

%PlateMapFile = fullfile(ExperimentFolder, ...
%    '20150727_MDAMB231_DoxorubicinCompartmentModel_PlateMap.csv');
%PlateMapFile = fullfile(ExperimentFolder, ...
%    '20150813_MDAMB231_DoxorubicinCompartmentModel_PlateMap.csv');
%PlateMapFile = fullfile(ExperimentFolder, ...
%    '20150902_MDAMB231_DoxorubicinCompartmentModel_PlateMap.csv');
%PlateMapFile = fullfile(ExperimentFolder, ...
%    '20151119_L12FluorescenceTest_PlateMap.csv');
PlateMapFile = fullfile(ExperimentFolder, ...
    '20160315_FixDoxSignalExperiment_PlateMap.csv');
%PlateMapFile = fullfile(ExperimentFolder, ...
%    '20160315_FixDoxSignalExperiment_002_PlateMap.csv');

% Read in plate map file to get concentrations and wells of interest
isDynamic = false;
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

c_foi = foi{21};
[allTimes_hrs, numImage_perFolder] = createTimeVector(DataFolders,c_foi);

%% Read and process images

% fluorIntenMatrix = number of wells x number of images x 4
% 4: mean background, mean cell, mean control, std background
% need to include control and cell wells when making fluorIntenMatrix
fluorIntenMatrix = zeros(8,sum(numImage_perFolder),4);

% processing
% subtract off corresponding control well from cell well
% Signal is summation of 2 compartments (cytoplasm, nuclear)

%for controlWellIter = 21:40%sum(controlWells)
for controlWellIter = 1:sum(controlWells)
    disp(controlWellIter)
    
    for dfIter = 1:length(DataFolders)
        
        controlFolder = fullfile(DataFolders{dfIter},controlFolders{controlWellIter});
        
        allOuttie =  zeros(numImage_perFolder(dfIter),4);
        
        for allIm = 0:numImage_perFolder(dfIter)-1
            
            flName = sprintf('Doxorubicin - n%06d.tif',allIm);
            control_flName = fullfile(controlFolder,flName);
            Im_control = double(imread(control_flName));
            
            outtie = [0 0 mean(Im_control(:)) std(Im_control(:))];
            allOuttie(allIm+1,:) = outtie;
            
        end
        
        startInsert = sum(numImage_perFolder(1:dfIter-1)) + 1;
        endInsert = sum(numImage_perFolder(1:dfIter));
        fluorIntenMatrix(controlWellIter, startInsert:endInsert, :) = allOuttie;
        
    end
    
end

%%
% look at fluorescent signal over time
%for aaa = 21:40

clrs = {'ro','go','bo','ko',...
    'r*','g*','b*','k*',...
    'rs','gs','bs','ks'};

figure(1022);clf;

uniqueConcs = unique(controlDrugConcs);
uniqueConcs = uniqueConcs(end:-1:1);

%for aaa = 1:30
for aaa = [1:4 10:14 20:24 30]
    figure(1022);
    hold on
    %if aaa<31
    %    plot(allTimes_hrs, fluorIntenMatrix(aaa,:,3),'bo')
    %else
    %    plot(allTimes_hrs, fluorIntenMatrix(aaa,:,3),'gs')
    %end
    
    ci = find(controlDrugConcs(aaa) == uniqueConcs);
    plot(allTimes_hrs, fluorIntenMatrix(aaa,:,3),clrs{ci})
    
    %plot(allTimes_hrs, fluorIntenMatrix(aaa,:,2),'ro')
    %plot(allTimes_hrs, fluorIntenMatrix(aaa,:,3),'ks')
    xlabel('Time (hrs)')
    ylabel('Intensity')
    axis([0 20 200 500])
    hold off
    
    %pause
end

%%

concMatrix = zeros(size(fluorIntenMatrix,1),size(fluorIntenMatrix,2),1);

figure(121);clf;
onlyThese = [1:4 10:14 20:24 30];
cdc = controlDrugConcs(onlyThese);

for tps = 1:size(fluorIntenMatrix,2)
    
    %cIM = fluorIntenMatrix(1:20,tps,3);
    %par = polyfit(controlDrugConcs(1:20),cIM,1);
    cIM = fluorIntenMatrix(onlyThese,tps,3);
    par = polyfit(cdc,cIM,1);
    
    figure(121);clf;
    hold on
    plot(cdc,cIM,'bo')
    %plot(controlDrugConcs(1:10),cIM(1:10),'gs')
    xdum = linspace(0,3000,10000);
    plot(xdum,par(1)*xdum+par(2),'r-')
    xlabel('Concentration (nM)')
    ylabel('Intensity')
    axis([-10 3000 220 400])
    
    concMatrix(:,tps) = (fluorIntenMatrix(:,tps,3)-par(2))./par(1);
    
    pause(.1)
end

concMatrix = concMatrix(1:30,:);
%concMatrix = concMatrix(1:20,:);

figure(12);clf;
imagesc(concMatrix)
xlabel('Time (hrs)')
ylabel('Plate Rows')