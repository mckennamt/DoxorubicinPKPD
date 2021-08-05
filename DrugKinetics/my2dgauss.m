

function outtie = my2dgauss(parms,X,Y,ImData,weightMatrix)

%%
%parms_est = [.5*size(ImData) 20 6 30 4000];
%parms = fitParms;
%parms(6) = 2800;

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
%keyboard
%figure(2);clf;
%subplot(121)
%imagesc(mulGau_2d)% figure(2);clf;
%subplot(121)
%imagesc(mulGau_2d)
%caxis([2000 4000])
%subplot(122)
%imagesc(ImData - mulGau_2d)
%caxis([-100 100])

outtie = (ImData(:) - mulGau_2d(:)).*weightMatrix(:);