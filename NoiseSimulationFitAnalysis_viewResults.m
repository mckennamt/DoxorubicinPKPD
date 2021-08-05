% full script to test fitting to simulated data with various noise levels

clear all;clc

load('NoiseSimulations.mat')


%%

% make a matrix to define average error rates

currFig = figure(1);clf

% nli = number of noise levels
allNoiseLevels = [0 .05 .1 .15 .2];

for nslvls = 1:nli
    

    holdTrue = allTrue{nslvls};
    holdPred = allPred{nslvls};
    %errorPred_percent = allError{nslvls};
    errorPred_percent = allCISize{nslvls};
    
    nrmlize = repmat(max(holdTrue), size(holdTrue,1), 1);
    %nrmlize = holdTrue;
    %nrmlize = (holdTrue + holdPred)./2
    errorPred_percent = 100* abs(holdTrue - holdPred)./nrmlize;
    
    
gridSize = 5;

kdRange = linspace(min(holdTrue(:,1)), max(holdTrue(:,1)), gridSize);
rRange = linspace(min(holdTrue(:,2)), max(holdTrue(:,2)), gridSize);

errorRates_kd = zeros(gridSize-1);
errorRates_r = zeros(gridSize-1);
for xi = 1:(gridSize-1)
    for yi = 1:(gridSize-1)
        
        ss_kp = holdTrue(:,1)>= kdRange(xi) & holdTrue(:,1) < kdRange(xi+1);
        ss_r = holdTrue(:,2)>= rRange(yi) & holdTrue(:,2) < rRange(yi+1);
        
        errorRates_kd(xi,yi) = mean(abs(errorPred_percent(ss_kp&ss_r, 1)));
        errorRates_r(xi,yi) = mean(abs(errorPred_percent(ss_kp&ss_r, 2)));
        
    end
end

%errorRates_kd(isnan(errorRates_kd)) = 1;


subplot(2,nli, nslvls);
[X,Y] = meshgrid(rRange(2:end),kdRange(2:end));
imagesc(rRange, kdRange, errorRates_kd)
%caxis([0 .1])
caxis([0 50])
xlabel('\it{r}');ylabel('\it{k_d_,_B}')
set(gca,'FontName','Times')
%title('kd error')

subplot(2,nli, nslvls+nli);
imagesc(rRange, kdRange, errorRates_r)
%title('r error')
%caxis([0 .1])
caxis([0 50])
xlabel('\it{r}');ylabel('\it{k_d_,_B}')
set(gca,'FontName','Times')
end

set(currFig, 'Position', [1550 320 1430 550])
colorbar