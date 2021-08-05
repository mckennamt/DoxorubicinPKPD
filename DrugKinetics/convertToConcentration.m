

function concMatrix = convertToConcentration(fluorIntenMatrix, allConcs, allTimes_hrs)

%keyboard
%%
concMatrix = zeros(size(fluorIntenMatrix,1),...
    size(fluorIntenMatrix,2), 2);

allConcs_input = allConcs(1:8);

for tps = 1:size(fluorIntenMatrix,2)
    
    cIM = fluorIntenMatrix(:,tps,3);
    
    if tps <= 4 %|| tps >=29
        allConcs_f = zeros(size(allConcs_input));
    else
        allConcs_f = allConcs_input;
    end
    
    if tps <= 15
par = polyfit(allConcs_f (cIM < 4000) ,cIM(cIM < 4000),1);
%par2 = polyfit(cIM(cIM < 4000),allConcs_f (cIM < 4000),1);
    end

% % cMod = linspace(0,10000);
% % fMod = par(1)*cMod + par(2);
% % figure(1);clf;
% % semilogx(allConcs_input, cIM,'r.')
% % hold on
% % semilogx(cMod,fMod,'b:')
% % %axis([0 11000 0 4500])
% % pause(1)

% figure(111);clf;
% hold on
% % %for aaa = 1:8
% % 
% % %end
% xDum = linspace(0,2600);
% plot(xDum,par(1)*xDum + par(2),'b-','LineWidth',2);
% errorbar(allConcs_input,cIM,fluorIntenMatrix(:,tps,4),'ro','MarkerSize',10,'MarkerFaceColor',[1 0 0],'LineWidth',2)
% ylabel('Intensity','FontSize',32)
% xlabel('Concentration (nM)','FontSize',32)
% axis([-100 2600 250 1000])
% set(gca,'FontSize',20)
% pause

concMatrix(:,tps,1) = (fluorIntenMatrix(:,tps,1)-par(2))./par(1);
concMatrix(:,tps,2) = (fluorIntenMatrix(:,tps,2)-par(2))./par(1);

end


% %%
% 
% figure(997);clf
% hold on
% eod = 7;
% for aaa = 3:8
% plot(allTimes_hrs(1:end-eod), concMatrix(aaa,1:end-eod,1),'g.-')
% plot(allTimes_hrs(1:end-eod), concMatrix(aaa,1:end-eod,2),'ro-')
% xlabel('Time (hrs)')
% ylabel('Intensity')
% axis([0 50 -100 10000])
% 
% pause
% end
% 
% hold off