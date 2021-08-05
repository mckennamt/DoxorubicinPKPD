



minVal = min(AIC_compartmentModel,[],2);
dAIC = AIC_compartmentModel - repmat(minVal, 1, size(AIC_compartmentModel,2));
expdAIC = exp(-0.5.*dAIC);
AIC_weights = expdAIC ./ repmat(sum(expdAIC,2),1,size(expdAIC,2))

figure(1010);clf;
mrks = {'ko','r*','bd'};
for iter = [3 1 2]
    plot(1:4, AIC_compartmentModel(:,iter),mrks{iter},'MarkerSize',10)
    %semilogy(1:4, dAIC(:,iter)+1,mrks{iter},'MarkerSize',12)
    hold on
end

axis([.5 4.5 -2000 500])

%legend('3 Compartment, 3 Paramter', '3 Compartment, 4 Parameter', '2 Compartment, 2 Parameter')
legend('Model a', 'Model b', 'Model c')

set(gca,'XTick',1:4,'XTickLabel',allCellLines,...
    'XTickLabelRotation',15,'FontSize',12)