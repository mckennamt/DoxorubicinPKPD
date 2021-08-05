
function [cellUptakeSignal_out, inputFcn_out] = correctUptakeSignal(cellUptakeSignal,...
    inputFcn, blankFcn, startDrug, c_endDrug)

% inputFcn = extracellular concentration
% daif = blank well with AIF
% cellUptakeSignal = Intracellular concentration
% blankFcn = well with no cells and no drug

%% modify intracellular curve

% remove background from image by subtracting off blankFcn both before and
% after drug is applied
cellUptakeSignal(1:startDrug) = cellUptakeSignal(1:startDrug) - blankFcn(1:startDrug);
cellUptakeSignal(c_endDrug+2:end) = cellUptakeSignal(c_endDrug+2:end) - blankFcn(c_endDrug+2:end);

% while drug is applied, subtract off extracellular signal
cellUptakeSignal(startDrug+1:c_endDrug+1) = cellUptakeSignal(startDrug+1:c_endDrug+1) - inputFcn(startDrug+1:c_endDrug+1);

% while drug is applied, subtract off extracellular signal
cellUptakeSignal(1:end) = cellUptakeSignal(:) - mean(cellUptakeSignal(1:startDrug));

% no drug prior to treatment
%cellUptakeSignal(1:startDrug) = cellUptakeSignal(1:startDrug) -...
%    mean(cellUptakeSignal(1:startDrug));
cellUptakeSignal(1:startDrug) = cellUptakeSignal(1:startDrug) -...
    mean(cellUptakeSignal(startDrug));

% concentration should always be >= 0
shiftUp = min(cellUptakeSignal(startDrug+1:c_endDrug+1));
cellUptakeSignal(startDrug+1:c_endDrug+1) = cellUptakeSignal(startDrug+1:c_endDrug+1)+abs(shiftUp).*(shiftUp<0);

%cellUptakeSignal = squeeze(concMatrix(wellIter,1:end,1))-squeeze(concMatrix(wellIter,1:end,2));
%cellUptakeSignal = squeeze(concMatrix(wellIter,1:end,1))-blankFcn;

%if min(cellUptakeSignal(c_endDrug+2:end)) <0
%    cellUptakeSignal(c_endDrug+2:end) = cellUptakeSignal(c_endDrug+2:end)+...
%        abs(min(inputFcn(c_endDrug+2:end)));
%end

%inputFcn(c_endDrug+2:end) = inputFcn(c_endDrug+2:end) + ...
%    abs(min(inputFcn(c_endDrug+2:c_endDrug+5)));

%% modify input function

% for the inputFcn (extracellular concentration) subtract off the no cell,
% no drug well
inputFcn = inputFcn-blankFcn;

% know 0 concentration before drug, shift down entire timecourse
inputFcn = inputFcn - mean(inputFcn(1:startDrug));

% know 0 concentration before drug, shift down entire timecourse
inputFcn(c_endDrug+2:end) = inputFcn(c_endDrug+2:end) - (inputFcn(c_endDrug+2));

%% if any value <0 shift up to zero

cellUptakeSignal(cellUptakeSignal<=0) = 0;
inputFcn(inputFcn<=0) = 0;

%% set output

cellUptakeSignal_out = cellUptakeSignal;
inputFcn_out = inputFcn;