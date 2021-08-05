

function [allTimes_hrs, numImage_perFolder] =...
    createTimeVector(DataFolders,c_foi)

numImage_perFolder = zeros(length(DataFolders),1);

for dfIter = 1:length(DataFolders)
    
    dirInfo = dir(fullfile(DataFolders{dfIter}, c_foi));
    dirNames = {dirInfo.name};
    doxFls = cellfun(@(x) ~isempty(strfind(x,'Doxorubicin - n')),dirNames,'UniformOutput',true);
    numImage_perFolder(dfIter) = sum(doxFls);
    
    dirDates = {dirInfo.date};
    doxDates = dirDates(doxFls);
    allTimes = cellfun(@(x) datenum(x,'dd-mmm-yyyy HH:MM:SS'),doxDates,'UniformOutput',true);
    
    if dfIter == 1
        firstExpTimepoint = allTimes(1);
        allTimes = allTimes - allTimes(1);
        allTimes_hrs = allTimes' * 24;
    else
        allTimes = allTimes - firstExpTimepoint;
        allTimes_hrs = [allTimes_hrs; allTimes' * 24];
    end
    
end