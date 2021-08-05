
function cellData = loadCellLineData(CellLine)

cellData.ExposureTimes = [6 12 24];

switch lower(CellLine)
    
    case 'mdamb231'
        
        % DrugAdded_TP = drug added immediately after
        %cellData.DrugAdded_Tp = 3;%20150923
        cellData.DrugAdded_Tp = 4;
        
        % Compartment model parameters for MDAMB231
        compartmentParms = [.06 .06/.13 .05];
        
        for expTime = 1:length(cellData.ExposureTimes)
            
            ExposureTime = cellData.ExposureTimes(expTime);
            
            if ExposureTime == 6
                
                %plateFile = 'Data/20150923_MDAMB231_Dox_6hr_PlateMap.csv';
                %dataFile = 'Data/20150923_MDAMB231_Dox_6hr_RecountData.xlsx';
                plateFile = 'Data/20160422_MDAMB231_Dox_6hr_PlateMap.csv';
                dataFile = 'Data/20160422_MDAMB231_Dox_6hr_Data.xlsx';
                isDynamic = true;
                
            elseif ExposureTime == 12
                
                %plateFile = 'Data/20150923_MDAMB231_Dox_12hr_PlateMap.csv';
                %dataFile = 'Data/20150923_MDAMB231_Dox_12hr_RecountData.xlsx';
                plateFile = 'Data/20160422_MDAMB231_Dox_12hr_PlateMap.csv';
                dataFile = 'Data/20160422_MDAMB231_Dox_12hr_Data.xlsx';
                isDynamic = true;
                
            elseif ExposureTime == 24
                
                %plateFile = 'Data/20150923_MDAMB231_Dox_24hr_PlateMap.csv';
                %dataFile = 'Data/20150923_MDAMB231_Dox_24hr_RecountData.xlsx';
                plateFile = 'Data/20160422_MDAMB231_Dox_24hr_PlateMap.csv';
                dataFile = 'Data/20160422_MDAMB231_Dox_24hr_Data.xlsx';
                isDynamic = true;
                
            end
            
            [plateMap,DrugTime] = importPlateMap(plateFile,isDynamic);
            cellData_c = importCellavistaData(dataFile,plateMap);
            cellData_c = calcDrugCurves(cellData_c, ExposureTime,...
                cellData.DrugAdded_Tp, compartmentParms);
            
            cellData = mergeCellData(cellData, cellData_c);
            
        end
        
        % observed some cell death in control wells, ignore that
        %cellData.goodControlData = 18;%20150923
        cellData.goodControlData = min(cellfun(@(x) length(x), cellData.times));
        
        %%For 20150923 dataset
        %%Incorrect objective used for timepoints 2 and 3
        %%Drug added after TP 5 (TP 3 when 2 and 3 are removed)
        %cellData.times = cellData.times([1 4:end]);
        %cellData.NumNuclei = cellData.NumNuclei(:,[1 4:end]);
        %cellData.NucDensity = cellData.NucDensity(:,[1 4:end]);
        
        
    case 'sum149'
        
        % DrugAdded_TP = drug added immediately after
        cellData.DrugAdded_Tp = 4;
        
        % Compartment model parameters for MDAMB231
        compartmentParms = [.06 .06/.13 .05];
        
        for expTime = 1:length(cellData.ExposureTimes)
            
            ExposureTime = cellData.ExposureTimes(expTime);
            
            if ExposureTime == 6
                
                plateFile = 'Data/20160127_SUM149_Dox_6hr_PlateMap.csv';
                %dataFile = 'Data/20160127_SUM149_Dox_6hr_Data.xlsx';
                dataFile = 'Data/20160127_SUM149_Dox_6hr_Data_Export0223.xlsx';
                isDynamic = false;
                
            elseif ExposureTime == 12
                
                plateFile = 'Data/20160127_SUM149_Dox_12hr_PlateMap.csv';
                %dataFile = 'Data/20160127_SUM149_Dox_12hr_Data.xlsx';
                dataFile = 'Data/20160127_SUM149_Dox_12hr_Data_Export0223.xlsx';
                isDynamic = false;
                
            elseif ExposureTime == 24
                
                plateFile = 'Data/20160127_SUM149_Dox_24hr_PlateMap.csv';
                %dataFile = 'Data/20160127_SUM149_Dox_24hr_Data.xlsx';
                dataFile = 'Data/20160127_SUM149_Dox_24hr_Data_Export0223.xlsx';
                isDynamic = false;
                
            end
            
            [plateMap,DrugTime] = importPlateMap(plateFile,isDynamic);
            cellData_c = importCellavistaData(dataFile,plateMap);
            cellData_c = calcDrugCurves(cellData_c, ExposureTime,...
                cellData.DrugAdded_Tp, compartmentParms);
            
            cellData = mergeCellData(cellData, cellData_c);
            
            
        end
        
        
        % observed some cell death in control wells, ignore that
        cellData.goodControlData = min(cellfun(@(x) length(x), cellData.times));
        
        
    case 'mdamb468'
        
        % DrugAdded_TP = drug added immediately after
        cellData.DrugAdded_Tp = 5;
        
        % Compartment model parameters for MDAMB231
        compartmentParms = [.0474 .0474/.4021 .0261];
        
        for expTime = 1:length(cellData.ExposureTimes)
            
            ExposureTime = cellData.ExposureTimes(expTime);
            
            if ExposureTime == 6
                
                plateFile = 'Data/20160221_MDAMB468_Dox_6hr_PlateMap.csv';
                dataFile = 'Data/20160221_MDAMB468_Dox_6hr_Data.xlsx';
                isDynamic = false;
                
            elseif ExposureTime == 12
                
                plateFile = 'Data/20160221_MDAMB468_Dox_12hr_PlateMap.csv';
                dataFile = 'Data/20160221_MDAMB468_Dox_12hr_Data.xlsx';
                isDynamic = false;
                
            elseif ExposureTime == 24
                
                plateFile = 'Data/20160221_MDAMB468_Dox_24hr_PlateMap.csv';
                dataFile = 'Data/20160221_MDAMB468_Dox_24hr_Data.xlsx';
                isDynamic = false;
                
            end
            
            [plateMap,DrugTime] = importPlateMap(plateFile,isDynamic);
            cellData_c = importCellavistaData(dataFile,plateMap);
            cellData_c = calcDrugCurves(cellData_c, ExposureTime,...
                cellData.DrugAdded_Tp, compartmentParms);
            
            cellData = mergeCellData(cellData, cellData_c);
            
        end
        
        % observed some cell death in control wells, ignore that
        cellData.goodControlData = min(cellfun(@(x) length(x), cellData.times))-1;
        
end



end

% combine all exposure time data for single cell line into one struct
function mergedCellData = mergeCellData(cellData, cellData_c)

mergedCellData = cellData;

cdf = fieldnames(cellData_c);

for iter = 1:length(cdf)
    currField = getfield(cellData_c,cdf{iter});
    
    if isfield(mergedCellData,cdf{iter})
        oldField = getfield(mergedCellData,cdf{iter});
        oldField{length(oldField)+1} = currField;
        mergedCellData = setfield(mergedCellData,cdf{iter},oldField);
        
        
    else
        mergedCellData = setfield(mergedCellData,cdf{iter},{currField});
    end
end
end

function cellData = calcDrugCurves(cellData, ExposureTime,...
    datp, compartmentParms)

uconcs = unique(cellData.drugConc);

% first column has Yd vector, second column has time vector corresponding
% to Yd vector
allYd_tv = cell(length(uconcs),3);
allYall = cell(length(uconcs),1);
dox_pss = zeros(length(uconcs),1);
auc = zeros(length(uconcs), 1);

for drugIter = 1:length(uconcs)
    
    wi = find(cellData.drugConc == uconcs(drugIter),1);
    
    % Pre-compute drug concentration curves
    Ts2 = cellData.times(datp:end);
    drugAdded = cellData.drugConc(wi);
    DrugSchedule.drugConc = [drugAdded 0];
    DrugSchedule.drugChangeTimes = [Ts2(1) Ts2(1)+ExposureTime];
    
    [daif,Yd, time_vector] = generateChemoAIF(Ts2(end), 1000, ...
        compartmentParms, DrugSchedule);
    daif(time_vector>Ts2(1)+ExposureTime) = 0;
    
    %%fixed to value after drug removed
    %eot = find(time_vector>=DrugSchedule.drugChangeTimes(2),1);
    %Yd(eot:end) = Yd(eot);
    
    allYd_tv{drugIter,1} = Yd;
    allYd_tv{drugIter,2} = time_vector-cellData.times(datp);
    allYd_tv{drugIter,3} = daif;
    
    %dox_pss(drugIter) = max(Yd);
    dox_pss(drugIter) = Yd(end);
    auc(drugIter) = drugAdded*ExposureTime;
    
    allYall{drugIter} = (cellData.NumNuclei(wi:(wi+5),1:end))';
end

Tspan = cellData.times-cellData.times(datp);

cellData.Tspan = Tspan;
cellData.allYall = allYall;
cellData.allYd_tv = allYd_tv;
cellData.dox_pss = dox_pss;
cellData.auc = auc;
cellData.uconcs = uconcs;


end