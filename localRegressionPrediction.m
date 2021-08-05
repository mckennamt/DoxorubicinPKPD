

for aaa = 5

figure(1);clf;

plot(log(currParms(:,1)), currParms(:,2).*currParms(:,3).^2,'k.')

 %f=fit([log(currParms(:,1)) ones(size(currParms(:,1)))],currParms(:,pIter+1), 'lowess')
 f=smooth(log(currParms(:,1)),currParms(:,2).*currParms(:,3).^2, aaa, 'loess');
 
 hold on
 plot(log(currParms(:,1)), f,'k-')
 
 plot(log(trueParms(:,1)),trueParms(:,2).*trueParms(:,3).^2,'ro')
 
 %pause(.5)
end

%%

figure(1);clf
subplot(121)
plot(log(trueParms(:,1)),trueParms(:,2).*trueParms(:,3).^2,'ro')
subplot(122)

plot(log(trueParms(:,1)),trueParms(:,3),'ro')

%%

%noisy = currParms(:,2).*currParms(:,3).^2;
noisy = currParms(:,2);

xy = [log(currParms(:,1)), noisy];
X = linspace(-2,7);
span = 100;

f = @(xy) mylowess(xy,X,span);
yboot2 = bootstrp(1000,f,xy)';

figure(1);clf;
meanloess = mean(yboot2,2);
h1 = line(X, meanloess,'color','k','linestyle','-','linewidth',1);

stdloess = std(yboot2,0,2);
h2 = line(X, meanloess+2*stdloess,'color','r','linestyle','--','linewidth',1);
h3 = line(X, meanloess-2*stdloess,'color','r','linestyle','--','linewidth',1);

hold 
%plot(log(trueParms(:,1)),trueParms(:,2).*trueParms(:,3).^2,'ro')
plot(log(trueParms(:,1)),trueParms(:,2),'bx')
plot(xy(:,1), xy(:,2),'ko')
