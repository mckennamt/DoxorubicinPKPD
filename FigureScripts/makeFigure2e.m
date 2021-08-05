
% Dynamic dox modeling
close all;clear all;clc

CellLine = 'SUM149';
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

analysisWells = [1:9 11:19];

%% Get time vectors and number of images

% Should be loaded with cell segmentation data
c_foi = foi{1};
[allTimes_hrs, numImage_perFolder] = createTimeVector(DataFolders,c_foi);
allTimes_hrs = allTimes_hrs - allTimes_hrs(startDrug+1);

%% for each timepoint, make concentration look up table

concFluorInten_Map = zeros(length(allTimes_hrs),2);

dfIter = 1;
allIm = 50;
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

bnds = tinv([0.025, 0.975],numel(Im_drug)*2-1)' * stdInten';%/sqrt(numel(Im_drug)*2);
ub = bnds(2,:);
lb = -bnds(1,:);

fs = 24;

cf = figure(1);clf;
%figure(121);clf
%plot(uConc,fluorInten,'r.')
ss = [1 6:length(uConc)];
hold on
cMod = linspace(0,2700);
fMod = linFit(1)*cMod + linFit(2);
errorbar(uConc(ss), fluorInten(ss), lb(ss), ub(ss),'ro','LineWidth',2)
plot(cMod,fMod,'k-','LineWidth',1.5);
axis([-100 2700 200 475]);
%%set(cf,'Position',[100 100 900 250])
%set(gca,'FontSize',18)
%pause(.5)
set(gca,'FontSize',fs)
set(cf,'Position',[80 1080 800 425])
cl = legend('Measurements', 'Linear Fit','Location','NorthWest');
set(cl,'FontSize',28)
