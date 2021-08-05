
keyboard
%%
parmSelect = 1;

mv = min(log(predictorVariables_test(:,1)));
Mv = max(log(predictorVariables_test(:,1)));
simlowess(:,1) = linspace(mv, Mv, 20);

mv = min(log(predictorVariables_test(:,2)));
Mv = max(log(predictorVariables_test(:,2)));
simlowess(:,2) = linspace(mv, Mv, 20);

smoothd_parms = mylowess([log(predictorVariables_training), allParms(:, parmSelect)],...
                    simlowess, validSpan);

cf = figure(121);clf;
semilogx(predictorVariables_training(:,1), allParms(:,parmSelect),'rx','MarkerSize',15)
hold on
plot(predictorVariables_test(:,1), predictedParms(:,parmSelect),'bo','MarkerSize',15)
plot(predictorVariables_test(:,1), predictedParms(:,parmSelect),'k-','LineWidth',2)

legend('Training Values',...
    'Predicted Values','Location','SouthEast')

set(cf,'Position',[100 100 800 500])
set(gca,'FontSize',14)
axis([.05 10^3.1 -.01 .03])