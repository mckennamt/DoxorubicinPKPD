
clear all;clc
cf = figure(1258);clf
set(cf,'Position',[1200 400 1000 200])

allCompart = {'SUM149CompartmentParameters.mat',...
    'MDAMB468CompartmentParameters.mat',...
    'MDAMB453CompartmentParameters.mat'};

for cellResultIter = 1:length(allCompart)

load(allCompart{cellResultIter})

allTimes_hrs_mod = allTimes_hrs + abs(min(allTimes_hrs));
options=optimset('TolFun',1e-6,'TolX',1e-6,'MaxIter',1e50,'Display','off','PrecondBandWidth',0);

% 3 compartment 3 parameter
% kpars =   [k1 k2 k3]
estim =     [2 0.05 .05];
lowBound =  [0.001  0.001 .001];
highBound = [10  5 5];
[parm,Rn,Res,whyExit,~,~,J] = lsqnonlin(@(kpars) runCompartmentModel_CellLine(kpars,concMatrix,...
    allTimes_hrs_mod, analysisWells, aifDuration,startDrug, endDrug),...
    estim, lowBound, highBound, options);

dof= (length(Res)-length(parm));
v = Rn/dof;
% the parameter covariance matrix
C = v*inv(J'*J);

% extract the standard errors for each fitted parameter
se_g = sqrt(diag(C));

% compute the 95% confidence intervals
tst = tinv(0.975,dof);
%disp(['95% confidence intervals of estimates of parameters'])
z_ci = [parm'-tst*se_g,parm',parm'+tst*se_g]



subplot(131)
hold on
p = 1;
%plot(1,z_ci(p,2,1),'ro')
errorbar(cellResultIter,z_ci(p,2,1),...
    z_ci(p,2,1)-z_ci(p,1,1),...
    z_ci(p,3,1)-z_ci(p,2,1),'ro')

subplot(132)
hold on
p = 2;
%plot(1,z_ci(p,2,1),'ko')
errorbar(cellResultIter,z_ci(p,2,1),...
    z_ci(p,2,1)-z_ci(p,1,1),...
    z_ci(p,3,1)-z_ci(p,2,1),'ko')

subplot(133)
hold on
p = 3;
%plot(1,z_ci(p,2,1),'bo')
errorbar(cellResultIter,z_ci(p,2,1),...
    z_ci(p,2,1)-z_ci(p,1,1),...
    z_ci(p,3,1)-z_ci(p,2,1),'bo')

end

subplot(131)
%axis([1e2 1e4 -.05 .2])
%axis([1e2 1e4 0 .2])

subplot(132)
%axis([1e2 1e4 -.1 1])
%axis([1e2 1e4 0 .1])

subplot(133)
%axis([1e2 1e4 0 .2])
%axis([1e2 1e4 0 2])