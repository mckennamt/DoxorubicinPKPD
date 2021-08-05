


cf = figure(1258);clf

set(cf,'Position',[1200 400 1000 200])

subplot(131)
semilogx(0,0)
hold on

subplot(132)
semilogx(0,0)
hold on

subplot(133)
semilogx(0,0)
hold on

uc = unique(allConcs(analysisWells));
z_conc = allConcs(analysisWells);


for iter = 1:length(uc)

    cc = uc(iter);

subplot(131)
semilogx(cc,mean(z_ci(1,2,z_conc==cc)),'ro')
errorbar(cc,mean(z_ci(1,2,z_conc==cc)),std(z_ci(1,2,z_conc==cc)),'ro')

subplot(132)
semilogx(cc,mean(z_ci(2,2,z_conc==cc)),'ko')
errorbar(cc,mean(z_ci(2,2,z_conc==cc)),std(z_ci(2,2,z_conc==cc)),'ko')

subplot(133)
semilogx(cc,mean(z_ci(3,2,z_conc==cc)),'bo')
errorbar(cc,mean(z_ci(3,2,z_conc==cc)),std(z_ci(3,2,z_conc==cc)),'bo')

end

subplot(131)
axis([1e2 1e4 -.05 .2])
%axis([1e2 1e4 0 .2])

subplot(132)
axis([1e2 1e4 -.1 1])
%axis([1e2 1e4 0 .1])

subplot(133)
%axis([1e2 1e4 0 .2])
%axis([1e2 1e4 0 2])



