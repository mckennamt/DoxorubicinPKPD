% normalize intensities


%%

for controlWellIter = 1
    disp(controlWellIter)
    
    for dfIter = 1:length(DataFolders)
        
        controlFolder = fullfile(DataFolders{dfIter},controlFolders{controlWellIter+8});
        cellFolder = fullfile(DataFolders{dfIter},cellFolders{controlWellIter});
        
        allImages = zeros(1536,2016,75);
        %allImages_cell = zeros(1536,2016,75);
        
        for allIm = 0:74%numImage_perFolder(dfIter)-1
            flName = sprintf('Doxorubicin - n%06d.tif',allIm);
            
            control_flName = fullfile(controlFolder,flName);
            cell_flName = fullfile(cellFolder,flName);
            allImages(:,:,allIm+1) = double(imread(control_flName));
            %allImages_cell(:,:,allIm+1) = double(imread(cell_flName));
            
            
        end
        
    end
end

%%

allImages_nc = double(allImages);
mv = max(allImages_nc(:));
allImages_nc = allImages_nc./mv;

v = 5;
MyTargets = [5:v:size(allImages,3)];


davec = 0:.01:1;

allImages_corrected = zeros(size(allImages_nc));
allImages_corrected(:,:,1:4) = allImages_nc(:,:,1:4);

for iter = 5:70
    
    iter
    
    %ti = find(MyTargets-iter>=0,1);
    %tn = MyTargets(ti);
    %wn = (v-(tn-iter))/v;
    
    %tnImage = allImages_nc(:,:,tn);
    %tnHist = hist(tnImage(:),davec);
    %Jn = histeq(allImages_nc(:,:,iter),tnHist);
    
    %tp = MyTargets(ti-1);
    %wp = (v-(iter-tp))/v;
    
    %tpImage = allImages_nc(:,:,tp);
    %tpHist = hist(tpImage(:),davec);
    %Jp = histeq(allImages_nc(:,:,iter),tpHist);
    
    validCount = 0;
    for maf = iter-3:iter+3
        %if (maf>=5 && maf<=28) && (iter>=5 && iter<=28)
        %    validCount = validCount + 1;
        %    tpImage = allImages_nc(:,:,maf);
        %    tpHist = hist(tpImage(:),davec);
        %    Jp = histeq(allImages_nc(:,:,iter),tpHist);
        %    allImages_corrected(:,:,iter) = Jp + ...
        %        allImages_corrected(:,:,iter);
        %elseif (maf<5 || maf>28) && (iter<5 || iter>28) && maf<76
        %    validCount = validCount + 1;
        %    tpImage = allImages_nc(:,:,maf);
        %    tpHist = hist(tpImage(:),davec);
        %    Jp = histeq(allImages_nc(:,:,iter),tpHist);
        %    allImages_corrected(:,:,iter) = Jp + ...
        %        allImages_corrected(:,:,iter);
        %end
        if (maf>=5) && (iter>=5)
            validCount = validCount + 1;
            tpImage = allImages_nc(:,:,maf);
            tpHist = hist(tpImage(:),davec);
            Jp = histeq(allImages_nc(:,:,iter),tpHist);
            allImages_corrected(:,:,iter) = Jp + ...
                allImages_corrected(:,:,iter);
        end
    end
    allImages_corrected(:,:,iter) = allImages_corrected(:,:,iter)./validCount;
    
    %allImages_corrected(:,:,iter) = wp*Jp + wn*Jn;
    
    %figure(12);clf;
    %subplot(121)
    %imagesc(allImages_corrected(:,:,iter)*mv)
    %caxis([0 1200])
    %title(num2str(iter))
    %subplot(122)
    %imagesc(allImages_cell(:,:,iter))
    %caxis([0 1200])
    %pause(.1)
    
    %subplot(131)
    %imagesc(allImages_nc(:,:,tp)*mv)
    %caxis([0 1200])
    %subplot(132)
    %imagesc(allImages_nc(:,:,tn)*mv)
    %caxis([0 1200])
    %subplot(133)
    %imagesc(allImages_corrected(:,:,iter)*mv)
    %caxis([0 1200])
    %pause(.1)
    
    
end
%%
theMeansCorrected = squeeze(mean(mean(allImages_corrected*mv,1),2));
theMeans = squeeze(mean(mean(allImages,1),2));
figure(100);clf;
hold on
plot(allTimes_hrs(1:70),theMeansCorrected(1:70),'bo-')
plot(allTimes_hrs(1:70),theMeans(1:70),'rs-')
legend('Corrected','Original','Location','Best')


% %%
%
%
% figure(1);clf;
% subplot(131)
% imagesc(allImages_nc(:,:,5)*mv)
% caxis([0 1200])
%
% subplot(132)
% imagesc(allImages_nc(:,:,6)*mv)
% caxis([0 1200])
%
% dubs = allImages_nc(:,:,6);
% [b,a] = hist(dubs(:),davec);
%
% dubs2 = allImages_nc(:,:,5);
% [b2,a2] = hist(dubs2(:),davec);
% %J = histeq(allImages(:,:,5)./mv,a);
% J = histeq(allImages_nc(:,:,5),b);
%
% figure(3);clf
% plot(a,b,'k.')
% hold on
% plot(a2,b2,'ro')
% [b3,a3] = hist(J(:),davec);
% plot(a3,b3,'gs')
%
% figure(1);
% subplot(133)
% imagesc(J*mv)
% caxis([0 1200])

