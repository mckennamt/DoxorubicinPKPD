
function [centroids,singleCell] = quickSeedGeneration(Im_cellTL,numCells, singleCell)

%numCells = 30;
%singleCell = Im_cellTL(168:205,187:228);
if isempty(singleCell)
h = figure(878);clf;
set(h,'name','Crop out single cell for cross-correlation','numbertitle','off') 
singleCell = imcrop(Im_cellTL./max(Im_cellTL(:)));
end

corrImage = 100*normxcorr2(singleCell,Im_cellTL);
cdn = imtophat(corrImage,strel('disk',10));
cdn = imdilate(cdn,ones(10,10));

bwd = imregionalmin(-cdn);
rp = regionprops(bwd,cdn,'MeanIntensity','Centroid');
mip = cat(1,rp.MeanIntensity);

centroids = cat(1,rp.Centroid);

pad = floor(size(singleCell)/2);
centroids(:,1) = centroids(:,1)-pad(1)*ones(size(centroids(:,1)));
centroids(:,2) = centroids(:,2)-pad(2)*ones(size(centroids(:,2)));

centroids = round(centroids);

imLimits = size(Im_cellTL);
removeMe = false(size(centroids,1),1);
removeMe(centroids(:,1)<1) = true;
removeMe(centroids(:,2)<1) = true;
removeMe(centroids(:,1)>imLimits(2)) = true;
removeMe(centroids(:,2)>imLimits(1)) = true;

centroids(removeMe,:) = [];
mip(removeMe) = [];

[~,si] = sort(mip,'descend');

autoCountCells = true;
while autoCountCells
    figure(878);
    imagesc(Im_cellTL)
    colormap gray
    hold on
    %keyboard
    curr_centroids = centroids(si(1:numCells),:);
    mySeed = zeros(size(Im_cellTL));
    seedInds = sub2ind(size(Im_cellTL),curr_centroids(:,2),curr_centroids(:,1));
    mySeed(seedInds) = 1;
    mySeed_dilated = imdilate(mySeed,ones(30,30));
    lab = label2rgb(mySeed_dilated,'jet');
    
    h = imshow(lab);
    set(h,'AlphaData',.5)
    
    waitforbuttonpress
    cb = double(get(gcf,'CurrentCharacter'));
    if cb == 31
        numCells = numCells-1;
    elseif cb == 30
        numCells = numCells+1;
    elseif cb == 13
        autoCountCells = false;
    end
    
    if numCells < 1
        numCells = 1;
    elseif numCells > size(centroids,1)
        numCells = size(centroids,1);
    end
    
end

centroids = centroids(si(1:numCells),:);
close
