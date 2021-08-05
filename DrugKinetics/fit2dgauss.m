
ImData = double(filtIm);

dubs = im2bw(ImData./max(ImData(:)),graythresh(ImData./max(ImData(:))));
weightMatrix = 1-double(dubs);

[X,Y] = meshgrid(1:1024,1:1024);

parms_est = [.5*size(ImData) 20 6 30 3800];

%ImData = double(Im_c);
lsqOpts = optimoptions('lsqnonlin','MaxFunEvals',200,'Display','iter');
[fitParms,rnorm,res] = ...
    lsqnonlin(@(p) my2dgauss(p,X,Y,ImData,weightMatrix),parms_est, ...
    [.25*size(ImData) 0 -Inf 0 0],...
    [.75*size(ImData) Inf Inf Inf Inf], lsqOpts);


%%
figure(2);clf;
subplot(121)
imagesc(ImData)
subplot(122)
parms = fitParms;
mean1=[parms(2); parms(1)];
cov1=[parms(3) parms(4);
    parms(4) parms(5)];

x1=[X(:) Y(:)]';
x1 = x1./100;
mean1 = mean1./100;


%multivar Gassiaan
mn=repmat(mean1,1,size(x1,2));

hldme = (x1-mn)'*inv(cov1);
aaa = sum(hldme.*(x1-mn)',2);
%mulGau= parms(6).*1/(2*pi*det(cov1)^(1/2))*...
%    exp(-0.5.*aaa);
mulGau= parms(6).*exp(-0.5.*aaa);


mulGau_2d = reshape(mulGau,size(X,1),size(X,2));

ou = histeq(ImData - mulGau_2d,512);

%imagesc(ImData - mulGau_2d)
imagesc(ou)
%imagesc(mulGau_2d)