% import plate map for cellavista data
% requires csv file

function [plateMap,DrugTime,cellType,folderList] = importPlateMap(fileName,isDynamic)

%keyboard
fid = fopen(fileName,'r');
headerLine = fgetl(fid);

% count number of columns
numComma = strfind(headerLine,',');
numColumns = length(numComma) + 1;

%load in all data
C = textscan(fid,['%s %s %s' repmat('%f',1,numColumns-3)],'Delimiter',',');
C{end} = [C{end};nan];
fclose(fid);

Column = cellfun(@(x) str2double(x(2:end)),C{1});
Row = cellfun(@(x) x(1),C{1});
folderList = C{1};
cellType = C{3};

DrugConc = C{7};

if isDynamic
    DrugTime = C{8};
else
    DrugTime = inf(size(DrugConc));
end

plateMap = {Row Column DrugConc DrugTime};