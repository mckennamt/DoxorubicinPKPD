

% function to load cellavista counts from ImageJ macro
function cellData = importImageJData(flname,plateMap,weightInfo)

[~,~,allCell] = xlsread(flname);
[~,c] = size(allCell);
emptyRow = sum(cell2mat(cellfun(@(x)sum(isnan(x))>0, ...
    allCell,'UniformOutput',false)),2) == c;
allCell = allCell(~emptyRow,:);

TimeStamps = logical(sum(cell2mat(cellfun(@(x)strcmp(x,'Timespan(h):'), ...
    allCell,'UniformOutput',false)),2));
times = [allCell{TimeStamps,2}]; % hours

measurementNumber = cell2mat(allCell([TimeStamps(2:end);false],1));

rowRange = find(TimeStamps,2);
allWells = cell2mat(allCell(rowRange(1)+1:rowRange(2)-2,1));
plateCol = unique(str2num(allWells(:,2:end)));
plateRow = cellstr(unique(allWells(:,1)));

NumNuclei = zeros(length(plateCol)*length(plateRow),max(measurementNumber));
drugConc = zeros(length(plateCol)*length(plateRow),1);
DrugTime = zeros(length(plateCol)*length(plateRow),1);

allDrug = plateMap{3};
allTimes = plateMap{4};

iter = 1;
for ci = 1:length(plateCol)
    for ri = 1:length(plateRow)
        c_col = plateCol(ci);
        c_row = plateRow{ri};
        
        luv = sprintf('%s%02d',c_row,c_col);
        cwell = cellfun(@(x) strcmp(x,luv), allCell(:,1));
        
        NumNuclei(iter,:) = [allCell{cwell,2}];%Estimated
        
        drugConc(iter) = allDrug(plateMap{2}==c_col & plateMap{1}==c_row);
        DrugTime(iter) = allTimes(plateMap{2}==c_col & plateMap{1}==c_row);
        
        iter = iter+1;
    end
end

cellData.times = times;
cellData.drugConc = drugConc;
cellData.NumNuclei = NumNuclei;
cellData.NumNuclei2 = [];
cellData.NucDensity = [];
cellData.DrugTime = DrugTime;

% load in last valid TP for each well
[~,~,allWeight] = xlsread(weightInfo.flname);
hoursColumns = cellfun(@(x) str2double(x(1:strfind(x, ' hour')-1)), allWeight(1,:));
cellData.LastValidTP = [allWeight{2:end,weightInfo.expTime == hoursColumns}]';


