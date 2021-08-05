
% Fit single set of compartment model parameters for each cell line

function R2 = runCompartmentModel_r2(parms,concMatrix,...
    allTimes_hrs_mod, analysisWells, aifDuration,startDrug,endDrug)

Y_meas = zeros(length(analysisWells), size(concMatrix,2));
Y_pred = zeros(size(Y_meas));
Y_weight = zeros(size(Y_meas));
cnt = 0;

for wellIter = analysisWells
    cnt = cnt+1;
    
    daif = squeeze(concMatrix(wellIter,1:end,3));
    inputFcn = squeeze(concMatrix(wellIter,1:end,2));
    blankFcn = squeeze(concMatrix(wellIter,1:end,4));
    
    if aifDuration(wellIter) == 6
        c_endDrug = endDrug(1);
    elseif aifDuration(wellIter) == 12
        c_endDrug = endDrug(2);
    end
    
    cellUptakeSignal = squeeze(concMatrix(wellIter,1:end,1));
    
    cellUptakeSignal(1:startDrug) = cellUptakeSignal(1:startDrug) - blankFcn(1:startDrug);
    cellUptakeSignal(c_endDrug+2:end) = cellUptakeSignal(c_endDrug+2:end) - blankFcn(c_endDrug+2:end);
    
    
    cellUptakeSignal(startDrug+1:c_endDrug+1) = cellUptakeSignal(startDrug+1:c_endDrug+1) - inputFcn(startDrug+1:c_endDrug+1);
    
    cellUptakeSignal(1:end) = cellUptakeSignal(:) - mean(cellUptakeSignal(1:startDrug));
    %cellUptakeSignal(startDrug+1:c_endDrug+1) = cellUptakeSignal(startDrug+1:c_endDrug+1)+abs(min(cellUptakeSignal));
    %cellUptakeSignal = squeeze(concMatrix(wellIter,1:end,1))-squeeze(concMatrix(wellIter,1:end,2));
    %cellUptakeSignal = squeeze(concMatrix(wellIter,1:end,1))-blankFcn;
    
    inputFcn = inputFcn-blankFcn;
    %inputFcn([1:4 20:end]) = inputFcn([1:4 20:end]) - daif([1:4 20:end]);
    %inputFcn(5:19) = inputFcn(5:19)-200;
    %inputFcn(c_endDrug+1:end) = 0;
    
    Y_meas(cnt,:) = cellUptakeSignal;
    %gtz = cellUptakeSignal + min(cellUptakeSignal);
    %Y_weight(cnt,:) = gtz./max(gtz);
    Y_weight(cnt,:) = 1;
    %Y_weight(cnt,c_endDrug+2:end) = 1;
    Y_pred(cnt,:) = modelCt3C3P(parms,inputFcn,allTimes_hrs_mod);
    
    %[~,Cout] = ode23(@(t,y) doxCarrierModel(t,y,parms,inputFcn, allTimes_hrs_mod),...
    %    allTimes_hrs_mod, [cellUptakeSignal(2) cellUptakeSignal(2)]);
    %Y_pred(cnt,:) = sum(Cout,2);
    
    %[~,Y_pred(cnt,:)] = ode45(@(t,y) twoPeakDoxModel(t,y,parms,inputFcn, allTimes_hrs_mod),...
    %    allTimes_hrs_mod, [cellUptakeSignal(2)]);
    
end

SS_Tot = Y_meas - repmat(mean(Y_meas,2),1,size(Y_meas,2));
SS_Tot = sum(SS_Tot(:).^2);

SS_Res = sum((Y_meas(:) - Y_pred(:)).^2);

R2 = 1 - SS_Res/SS_Tot;

