

function modSeeds = seedShift(prevImage,currentImage, currentSeeds)

prevIm_cropped = prevImage(100:end-100,100:end-100);
corrImage = 100*normxcorr2(prevIm_cropped,currentImage);
cdn = imtophat(corrImage,strel('disk',20));

% find maximum value of correlation
[~,mInd] = max(cdn(:));

% convert to coordinates, adapt for larger size of correlation image
[cx, cy] = ind2sub(size(corrImage),mInd);
pad = floor(size(prevIm_cropped)/2);
true_cx = cx-pad(1);
true_cy = cy-pad(2);

shiftVector = round([true_cx true_cy]) - size(currentImage)./2;

if sum(shiftVector > 20)==0
    
    modSeeds = currentSeeds + repmat(shiftVector([2 1]), size(currentSeeds,1),1);
    
    imLimits = size(currentImage);
    modSeeds(modSeeds(:,1)>imLimits(2) | modSeeds(:,1)<1,:) = [];
    modSeeds(modSeeds(:,2)>imLimits(1) | modSeeds(:,2)<1,:) = [];
else
    modSeeds = currentSeeds;
end

% figure(121)
% imagesc(cdn)
%
% figure(12122);clf;
% imagesc(currentImage)
% colormap gray
% prevSeeds = zeros(size(prevImage));
% %seedInds = sub2ind(size(prevImage),floor(modSeeds(:,2)),floor(modSeeds(:,1)));
% seedInds = sub2ind(size(prevImage),floor(currentSeeds(:,2)),floor(currentSeeds(:,1)));
% prevSeeds(seedInds) = 1;
% prevSeeds_dilated = imdilate(prevSeeds,ones(30,30));
%
% prevSeeds2 = zeros(size(prevImage));
% seedInds = sub2ind(size(prevImage),floor(modSeeds(:,2)),floor(modSeeds(:,1)));
% prevSeeds2(seedInds) = 1;
% prevSeeds2_dilated = imdilate(prevSeeds2,ones(30,30));
%
% hold on
% h = imshow(cat(3,ones(size(prevImage)),zeros(size(prevImage)),zeros(size(prevImage))));
% set(h,'AlphaData',.2.*prevSeeds_dilated);
%
% h2 = imshow(cat(3,zeros(size(prevImage)),zeros(size(prevImage)),ones(size(prevImage))));
% set(h2,'AlphaData',.2.*prevSeeds2_dilated);