%close all
%figure(1);clf;

%allExposureTimes = [6 12 24];
%allCellLine = {'MDAMB231','SUM149','MDAMB468'};
%cellLineSelect = 2;
%CellLine = allCellLine{cellLineSelect};


for cellLineSelect = 1:3
    figure(cellLineSelect);clf;
    set(gcf,'name',allCellLine{cellLineSelect},'numbertitle','off') 
    
currentP = collect_Fits{cellLineSelect};

daColors = {'ro-','gs-','b*-'};
extc = 0;
for ExposureTime = 1:3
    extc = extc+1;
    for functionSelectIter = 3
        
        holdParms = currentP{ExposureTime, functionSelectIter};
        
        %%first concentration is untreated, useless parameter
        %fc = 2;
        
        % AUC or CB (bound drug) or just concentration (no time considered)?
        %uc = cellfun(@(x) max(x), allYd_tv(fc:end,1));%CB
        %uc = cellfun(@(x) max(x), allYd_tv(fc:end,3))*ExposureTime;%AUC
        %uc = uconcs(fc:end);% just concentration
        uc = allCellData{cellLineSelect}.dox_pss{ExposureTime}(2:end);
        
        
        for iter = 1:size(holdParms,2)
        cp = holdParms(:,iter);
        
        %lb = squeeze(allParms_ci(iter,1,2:end));
        %ub = squeeze(allParms_ci(iter,2,2:end));
        
        subplot(1,size(holdParms,2),iter)
        loglog(uc,cp,daColors{extc},'MarkerSize',2,'LineWidth',2)
        set(gca,'Xlim',[.1 400])
        hold on
        %errorbar(uc,cp,abs(cp-lb),abs(ub-cp),daColors{extc})
        
        end
        
        
    end
end

legend('6 hour','12 hour','24 hour')
%axis([40 10^5 -5 60])
%axis([40 10^5 -.05 .2])
set(gca,'FontSize',12)

end