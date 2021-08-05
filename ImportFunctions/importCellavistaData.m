

% function to load cellavista data
function cellData = importCellavistaData(flname,plateMap)

[~,~,allCell] = xlsread(flname);
[~,c] = size(allCell);
emptyRow = sum(cell2mat(cellfun(@(x)sum(isnan(x))>0, ...
    allCell,'UniformOutput',false)),2) == c;
allCell = allCell(~emptyRow,:);

TimeStamps = logical(sum(cell2mat(cellfun(@(x)strcmp(x,'Timespan (h):'), ...
    allCell,'UniformOutput',false)),2));
times = [allCell{TimeStamps,2}]; % hours

measurementNumber = cell2mat(allCell([TimeStamps(2:end);false],1));

rowRange = find(TimeStamps,2);
plateCol = unique(cell2mat(allCell(rowRange(1)+1:rowRange(2)-2,1)));
plateRow = unique(allCell(rowRange(1)+1:rowRange(2)-2,2));

NucDensity = zeros(length(plateCol)*length(plateRow),max(measurementNumber));
NumNuclei = zeros(length(plateCol)*length(plateRow),max(measurementNumber));
NumNuclei2 = zeros(length(plateCol)*length(plateRow),max(measurementNumber));
drugConc = zeros(length(plateCol)*length(plateRow),1);
DrugTime = zeros(length(plateCol)*length(plateRow),1);

allDrug = plateMap{3};
allTimes = plateMap{4};
%keyboard
iter = 1;
for ci = 1:length(plateCol)
    for ri = 1:length(plateRow)
        c_col = plateCol(ci);
        c_row = plateRow{ri};
        
        cs = cell2mat(cellfun(@(x) sum(x==c_col),allCell(:,1),'UniformOutput',false));
        cr = cell2mat(cellfun(@(x) sum(x==c_row),allCell(:,2),'UniformOutput',false));        
        
        NucDensity(iter,:) = [allCell{cs&cr,3}];
        NumNuclei(iter,:) = [allCell{cs&cr,5}];%Estimated
        NumNuclei2(iter,:) = [allCell{cs&cr,7}];%Actually counted
        
        drugConc(iter) = allDrug(plateMap{2}==c_col & plateMap{1}==c_row);
        DrugTime(iter) = allTimes(plateMap{2}==c_col & plateMap{1}==c_row);
        
        iter = iter+1;
    end 
end

cellData.times = times;
cellData.drugConc = drugConc;
cellData.NumNuclei = NumNuclei;
cellData.NumNuclei2 = NumNuclei2;
cellData.NucDensity = NucDensity;
cellData.DrugTime = DrugTime;