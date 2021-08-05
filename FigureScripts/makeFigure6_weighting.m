
keyboard
%%

minVal = min(collectAICs,[],2);
dAIC = collectAICs - repmat(minVal, 1, size(collectAICs,2));
expdAIC = exp(-0.5.*dAIC);
AIC_weights = expdAIC ./ repmat(sum(expdAIC,2),1,size(expdAIC,2));

% fit logistic model to these weights
b_coeff = glmfit(trainVar, [AIC_weights(:,2) ones(size(AIC_weights(:,2)))],'binomial');

% calculate predicted value for logistic model with coefficients b_coeff at
% the points defined in predVar (for the predicted model)
pred_weights = glmval(b_coeff,predVar,'logit');

cf = figure(11000);clf;

semilogx(exp(trainVar(:,1)), AIC_weights(:,2),'rx','MarkerSize',15)
hold on
plot(exp(predVar(:,1)), pred_weights,'bo','MarkerSize',15)
plot(exp(predVar(:,1)), pred_weights,'k-','Linewidth',2)
%plot(predVar(:,1), 1-glmval(b_coeff,predVar,'logit'),'bo--','Linewidth',2)

%axis([.1 1e3 -.01 1.05])

legend('Training Weights',...
    'Predicted Weights','Location','SouthWest')

set(cf,'Position',[100 100 800 500])
set(gca,'FontSize',14)
axis([.05 10^3.1 -.05 1.05])