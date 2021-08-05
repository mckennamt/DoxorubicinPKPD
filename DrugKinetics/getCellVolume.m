

function averageVolume = getCellVolume(cellLine)

% calculate cell volume (um^3)
% assuming prolate spheroid (ellipse rotated about major axis)

cellSizeFile = [cellLine '_CellSizeMeasurements.xlsx'];
% measurements in micrometers (um)

[xl,xl_text] = xlsread(cellSizeFile);
majorAxCol = strcmp('Major', xl_text(1,:));
minorAxCol = strcmp('Minor', xl_text(1,:));
majorAxis = xl(:,majorAxCol);
minorAxis = xl(:,minorAxCol);

%%ellipsoid model
%volumeEstimate = 4/3 * pi * (minorAxis./2).^2 .* majorAxis./2;

%fixed height model
volumeEstimate = pi .* (minorAxis./2.* majorAxis./2) .* 1;% fixed height


averageVolume = mean(volumeEstimate);