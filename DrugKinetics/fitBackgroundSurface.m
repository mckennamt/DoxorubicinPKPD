

function minThis = fitBackgroundSurface(parms, ImData)

centerIm = zeros(size(ImData));
parms = floor(parms);
ci = sub2ind(size(ImData),parms(1),parms(2));
centerIm(ci) = 1;
r = double(bwdist(centerIm));

% define model
% parabola
C = [r(:).^2 r(:).^1 ones(size(r(:)))];
parms = lsqlin(C,double(ImData(:)),[],[]);
outtie = parms(1) * r.^2 + ...
    parms(2) * r.^1 + ...
    parms(3);

figure(112)
subplot(121)
imagesc(outtie)
subplot(122)
imagesc(double(-log(double(ImData)./max(ImData(:))))./outtie);caxis([-10 10])

minThis = double(-log(double(ImData)./max(ImData(:)))) - outtie;
minThis = minThis(:);


% % %correct background on fluorescent images!
% % 
% % % I = I_true * T
% % % I/T * I_median = I_true
% % 
% % 
% % ImData = double(Im_c);
% % lsqOpts = optimoptions('lsqnonlin','MaxFunEvals',300,'Display','iter');
% % [fitParms,rnorm,res] = ...
% %     lsqnonlin(@(p) fitBackgroundSurface(p,ImData),.5*size(ImData), ...
% %     floor(.25*size(ImData)), ceil(.75*size(ImData)), lsqOpts);
% % 
% % fitParms = floor(fitParms);
% % centerIm = zeros(size(ImData));
% % ci = sub2ind(size(ImData),fitParms(1),fitParms(2));
% % centerIm(ci) = 1;
% % r = double(bwdist(centerIm));
% % 
% % % define model
% % C = [r(:).^2 r(:).^1 ones(size(r(:)))];
% % parms = lsqlin(C,double(ImData(:)) - ones(size(r(:))),[],[]);
% % 
% % outtie = parms(1) * r.^2 + ...
% %     parms(2) * r.^1 + ...
% %     parms(3);
% % 
% % [X,Y] = meshgrid(1:size(ImData,2),1:size(ImData,1));
% % aaa = fit([X(:),Y(:)],double(ImData(:)), 'gauss2')
% % 
% % 
% % figure(112);clf;
% % subplot(121)
% % imagesc(outtie)
% % subplot(122)
% % imagesc(double(ImData))
% % 
% % figure(11)
% % subplot(133)
% % imagesc(double(Im_c) - outtie);caxis([0 10])


