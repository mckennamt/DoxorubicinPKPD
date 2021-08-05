
parentDataFolder = '/media/matthew/DataDisk/';

switch lower(CellLine)
    
    case 'sum149'
        noCurrData = exist('currDate');
        if ~noCurrData
            
            % SUM149; 3/17/2016
            ExperimentFolder = {parentDataFolder 'DrugKinetics_Data/20160317_DoxCompartmentModel_SUM149'};
            DataFolders{1} =  fullfile(ExperimentFolder, '2016-03-17_000');
            PlateMapFile = fullfile(ExperimentFolder, ...
                '20160317_SUM149_DoxorubicinCompartmentModel_PlateMap.csv');
            startDrug = 4;%Image file number (first image with drug on cells; zero-index)
            endDrug = [18 32];%Image file number (last image with drug on cells; zero-index)
            
        elseif currDate == 20161018
            
            % SUM149; 10/18/2016
            ExperimentFolder = {parentDataFolder 'DrugKinetics_Data/20161018_DoxCompartmentModel_SUM149'};
            DataFolders{1} =  fullfile(ExperimentFolder, '2016-10-18_000');
            PlateMapFile = fullfile(ExperimentFolder, ...
                '20161018_SUM149_DoxorubicinCompartmentModel_PlateMap.csv');
            startDrug = 4;%Image file number (first image with drug on cells; zero-index)
            endDrug = [22 39];%Image file number (last image with drug on cells; zero-index)
            
        elseif currDate == 20161022
            
            % SUM149; 10/22/2016
            ExperimentFolder = {parentDataFolder 'DrugKinetics_Data/20161022_DoxCompartmentModel_SUM149'};
            DataFolders{1} =  fullfile(ExperimentFolder, '2016-10-22_000');
            PlateMapFile = fullfile(ExperimentFolder, ...
                '20161022_SUM149_DoxorubicinCompartmentModel_PlateMap.csv');
            startDrug = 4;%Image file number (first image with drug on cells; zero-index)
            endDrug = [21 38];%Image file number (last image with drug on cells; zero-index)
        end
        
    case 'mdamb468'
        % MDAMB468; 3/18/2016
        ExperimentFolder = {parentDataFolder 'DrugKinetics_Data/20160318_DoxCompartmentModel_MDAMB468'};
        DataFolders{1} =  fullfile(ExperimentFolder, '2016-03-18_000');
        PlateMapFile = fullfile(ExperimentFolder, ...
            '20160318_MDAMB468_DoxorubicinCompartmentModel_PlateMap.csv');
        startDrug = 6;%Image file number (first image with drug on cells; zero-index)
        endDrug = [20 34];%Image file number (last image with drug on cells; zero-index)
        
    case 'mdamb453'
        
        % MDAMB453; 3/20/2016
        ExperimentFolder = {parentDataFolder 'DrugKinetics_Data/20160320_DoxCompartmentModel_MDAMB453'};
        DataFolders{1} =  fullfile(ExperimentFolder, '2016-03-20_000');
        PlateMapFile = fullfile(ExperimentFolder, ...
            '20160320_MDAMB453_DoxorubicinCompartmentModel_PlateMap.csv');
        startDrug = 5;%Image file number (first image with drug on cells; zero-index)
        endDrug = [19 33];%Image file number (last image with drug on cells; zero-index)
        
    case 'mdamb231'
        
        %%MDAMB231; 4/12/2016
        %ExperimentFolder = {parentDataFolder 'DrugKinetics_Data/20160412_DoxCompartmentModel_MDAMB231'};
        %DataFolders{1} =  fullfile(ExperimentFolder, '2016-04-12_000');
        %PlateMapFile = fullfile(ExperimentFolder, ...
        %    '20160412_MDAMB231_DoxorubicinCompartmentModel_PlateMap.csv');
        %startDrug = 4;%Image file number (first image with drug on cells; zero-index)
        %endDrug = [18 32];%Image file number (last image with drug on cells; zero-index)
        
        % MDAMB231; 9/8/2016
        ExperimentFolder = {parentDataFolder 'DrugKinetics_Data/20160908_DoxCompartmentModel_MDAMB231'};
        DataFolders{1} =  fullfile(ExperimentFolder, '2016-09-08_000');
        PlateMapFile = fullfile(ExperimentFolder, ...
            '20160908_MDAMB231_DoxorubicinCompartmentModel_PlateMap.csv');
        startDrug = 5;%Image file number (first image with drug on cells; zero-index)
        endDrug = [22 40];%Image file number (last image with drug on cells; zero-index)
        
end