
function cellData = loadCellLineData_IJ(CellLine)

cellData.ExposureTimes = [6 12 24];

switch lower(CellLine)
    
    case 'mdamb231'
        
        % DrugAdded_TP = drug added immediately after
        %cellData.DrugAdded_Tp = 3;%20150923
        cellData.DrugAdded_Tp = 4;
        
        % Compartment model parameters for MDAMB231
        compartmentParms = [.1867 .1867/.2410 .0865];
        
        for expTime = 1:length(cellData.ExposureTimes)
            
            ExposureTime = cellData.ExposureTimes(expTime);
            
            plateFile = ['Data/20160422_MDAMB231_Dox_' num2str(ExposureTime) 'hr_PlateMap.csv'];
            dataFile = ['Data/CVProcessing_ImageJ/20160422_MDAMB231_Dox_' num2str(ExposureTime) 'hr.xlsx'];
            weightInfo.flname = 'Data/CVProcessing_ImageJ/MDAMB231_CountWeights.xlsx';
            weightInfo.expTime = ExposureTime;
            isDynamic = false;
            
            [plateMap,DrugTime] = importPlateMap(plateFile,isDynamic);
            %cellData_c = importCellavistaData(dataFile,plateMap);
            cellData_c = importImageJData(dataFile,plateMap,weightInfo);
            cellData_c = calcDrugCurves(cellData_c, ExposureTime,...
                cellData.DrugAdded_Tp, compartmentParms);
            
            cellData = mergeCellData(cellData, cellData_c);
            
        end
        
        % observed some cell death in control wells, ignore that
        %cellData.goodControlData = 18;%20150923
        cellData.goodControlData = min(cellfun(@(x) length(x), cellData.times));
        % replaced with cellData.LastValidTP
        
        
    case 'sum149'
        
        % DrugAdded_TP = drug added immediately after
        cellData.DrugAdded_Tp = 4;
        
        % Compartment model parameters for SUM149
        compartmentParms = [.1812 .1812/.1434 .2095];
        
        for expTime = 1:length(cellData.ExposureTimes)
            
            ExposureTime = cellData.ExposureTimes(expTime);
            
            plateFile = ['Data/20160127_SUM149_Dox_' num2str(ExposureTime) 'hr_PlateMap.csv'];
            dataFile = ['Data/CVProcessing_ImageJ/20160127_SUM149_Dox_' num2str(ExposureTime) 'hr.xlsx'];
            weightInfo.flname = 'Data/CVProcessing_ImageJ/SUM149_CountWeights.xlsx';
            weightInfo.expTime = ExposureTime;
            isDynamic = false;
            
            [plateMap,DrugTime] = importPlateMap(plateFile,isDynamic);
            %cellData_c = importCellavistaData(dataFile,plateMap);
            cellData_c = importImageJData(dataFile,plateMap,weightInfo);
            cellData_c = calcDrugCurves(cellData_c, ExposureTime,...
                cellData.DrugAdded_Tp, compartmentParms);
            
            cellData = mergeCellData(cellData, cellData_c);
            
            
        end
        
        
        % observed some cell death in control wells, ignore that
        cellData.goodControlData = min(cellfun(@(x) length(x), cellData.times));
        % replaced with cellData.LastValidTP
        
        
    case 'mdamb468'
        
        % DrugAdded_TP = drug added immediately after
        cellData.DrugAdded_Tp = 5;
        
        % Compartment model parameters for MDAMB468
        compartmentParms = [.0534 .0534/.2611 .0728];
        
        for expTime = 1:length(cellData.ExposureTimes)
            
            ExposureTime = cellData.ExposureTimes(expTime);
            
            plateFile = ['Data/20160221_MDAMB468_Dox_' num2str(ExposureTime) 'hr_PlateMap.csv'];
            dataFile = ['Data/CVProcessing_ImageJ/20160221_MDAMB468_Dox_' num2str(ExposureTime) 'hr.xlsx'];
            weightInfo.flname = 'Data/CVProcessing_ImageJ/MDAMB468_CountWeights.xlsx';
            weightInfo.expTime = ExposureTime;
            isDynamic = false;
            
            [plateMap,DrugTime] = importPlateMap(plateFile,isDynamic);
            %cellData_c = importCellavistaData(dataFile,plateMap);
            cellData_c = importImageJData(dataFile,plateMap,weightInfo);
            cellData_c = calcDrugCurves(cellData_c, ExposureTime,...
                cellData.DrugAdded_Tp, compartmentParms);
            
            cellData = mergeCellData(cellData, cellData_c);
            
        end
        
        % observed some cell death in control wells, ignore that
        cellData.goodControlData = min(cellfun(@(x) length(x), cellData.times))-1;
        % replaced with cellData.LastValidTP
        
    case 'mdamb453'
        
        % DrugAdded_TP = drug added immediately after
        cellData.DrugAdded_Tp = 5;
        
        % Compartment model parameters for MDAMB453
        compartmentParms = [.0390 .0390/.4353 .0806];
        
        for expTime = 1:length(cellData.ExposureTimes)
            
            ExposureTime = cellData.ExposureTimes(expTime);
            
            plateFile = ['Data/20160313_MDAMB453_Dox_' num2str(ExposureTime) 'hr_PlateMap.csv'];
            dataFile = ['Data/CVProcessing_ImageJ/20160313_MDAMB453_Dox_' num2str(ExposureTime) 'hr.xlsx'];
            weightInfo.flname = 'Data/CVProcessing_ImageJ/MDAMB453_CountWeights.xlsx';
            weightInfo.expTime = ExposureTime;
            isDynamic = false;
            
            [plateMap,DrugTime] = importPlateMap(plateFile,isDynamic);
            %cellData_c = importCellavistaData(dataFile,plateMap);
            cellData_c = importImageJData(dataFile,plateMap,weightInfo);
            cellData_c = calcDrugCurves(cellData_c, ExposureTime,...
                cellData.DrugAdded_Tp, compartmentParms);
            
            cellData = mergeCellData(cellData, cellData_c);
            
        end
        
        % observed some cell death in control wells, ignore that
        cellData.goodControlData = min(cellfun(@(x) length(x), cellData.times))-1;
        % replaced with cellData.LastValidTP
        
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
allYweights = cell(length(uconcs),1);
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
    
    if isfield(cellData,'LastValidTP')
        allYweights{drugIter} = (cellData.LastValidTP(wi:(wi+5)))';
    else
        allYweights{drugIter} = [];
    end
end

Tspan = cellData.times-cellData.times(datp);

cellData.Tspan = Tspan;
cellData.allYall = allYall;
cellData.allYweights = allYweights;
cellData.allYd_tv = allYd_tv;
cellData.dox_pss = dox_pss;
cellData.auc = auc;
cellData.uconcs = uconcs;


end