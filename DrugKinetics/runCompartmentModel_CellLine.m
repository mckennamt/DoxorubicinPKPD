
% Fit single set of compartment model parameters for each cell line

function [Y,R2] = runCompartmentModel_CellLine(parms,concMatrix,...
    allTimes_hrs_mod, analysisWells, aifDuration,startDrug,endDrug)

Y_meas = zeros(length(analysisWells), size(concMatrix,2));
Y_pred = zeros(size(Y_meas));
Y_weight = zeros(size(Y_meas));
cnt = 0;

for wellIter = analysisWells
    cnt = cnt+1;
    
    cellUptakeSignal = squeeze(concMatrix(wellIter,1:end,1));
    inputFcn = squeeze(concMatrix(wellIter,1:end,2));
    daif = squeeze(concMatrix(wellIter,1:end,3));
    blankFcn = squeeze(concMatrix(wellIter,1:end,4));
    
    if aifDuration(wellIter) == 6
        c_endDrug = endDrug(1);
    elseif aifDuration(wellIter) == 12
        c_endDrug = endDrug(2);
    else
        c_endDrug = endDrug(1);
    end
    
    % clean up uptake curves
    [cellUptakeSignal, inputFcn] = correctUptakeSignal(cellUptakeSignal,...
        inputFcn, blankFcn, startDrug, c_endDrug);
    
    Y_meas(cnt,:) = cellUptakeSignal;
    
    % weighting data
    AUC=trapz(allTimes_hrs_mod,inputFcn);
    Y_weight(cnt,:) = 1*(cellUptakeSignal>=0);% & inputFcn>0);
    Y_weight(cnt,:)=Y_weight(cnt,:)/AUC;
    
    % select model to use
    Y_pred(cnt,:) = modelCt3C3P(parms,inputFcn,allTimes_hrs_mod);
    %Y_pred(cnt,:) = modelCt3C4P(parms,inputFcn,allTimes_hrs_mod);
    %Y_pred(cnt,:) = modelCt2C2P(parms,inputFcn,allTimes_hrs_mod);
    
    % hand-coded 3 compartment model (equivalent to modelCt3C3P)
    %[~, Cout] = ode45(@(t,y) chemoCompartmentModel_odes(t,y,...
    %    parms,inputFcn,allTimes_hrs_mod),...
    %    allTimes_hrs_mod,[cellUptakeSignal(2) cellUptakeSignal(2)]);
    %Y_pred(cnt,:) = sum(Cout,2);
    
    %[~,Cout] = ode45(@(t,y) doxCarrierModel(t,y,parms,...
    %    inputFcn, allTimes_hrs_mod),...
    %    allTimes_hrs_mod, [0 0]);
    %Y_pred(cnt,:) = sum(Cout,2);
    
    %[~,Y_pred(cnt,:)] = ode45(@(t,y) twoPeakDoxModel(t,y,parms,inputFcn, allTimes_hrs_mod),...
    %    allTimes_hrs_mod, [cellUptakeSignal(2)]);
    
end

% calculate r^2 of fit
Y_weight_mod = Y_weight;
Y_weight_mod(Y_weight_mod>0) = 1;
SS_Tot = Y_meas - repmat(mean(Y_meas,2),1,size(Y_meas,2));
SS_Tot = sum(SS_Tot(:).^2.*Y_weight_mod(:));
SS_Res = sum((Y_meas(:) - Y_pred(:)).^2.*Y_weight_mod(:));
R2 = 1 - SS_Res/SS_Tot;

% calculate mean percent error
%Y_weight_mod
R2 = mean(abs((Y_meas(Y_meas>0) - Y_pred(Y_meas>0))./(Y_meas(Y_meas>0)).*(Y_weight_mod(Y_meas>0))));
%R2 = mean(abs((Y_meas(Y_weight_mod>0) - Y_pred(Y_weight_mod>0))./(Y_meas(Y_weight_mod>0)).*(Y_weight_mod(Y_weight_mod>0))));

Y_weight(:,end) = 0;
Y_pred = Y_pred(:);
Y_meas = Y_meas(:);
Y_weight = Y_weight(:);
%Y = ((Y_meas - Y_pred)./(abs(Y_meas)+1)).*Y_weight;

% Are there alternative weighting methods to account for the noise in the
% system (particularly at the lower concentrations with minimal drug
% uptake)
Y = (Y_meas - Y_pred).*sqrt(Y_weight);