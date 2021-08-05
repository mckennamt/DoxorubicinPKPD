
function [seg_im,currentSeeds,proceedToTL] = segment_H2BRFPImages(Im_cellH2B)

% get initial seeds from simple segmentation information
proceedToTL = false;

%% rolling ball filter to subtract out background
diskSize = 20;
rbb = imfilter(Im_cellH2B, fspecial('disk', diskSize), 'replicate','same');

fluorIm = Im_cellH2B - rbb;
fluorIm = fluorIm./max(fluorIm(:));

%% threshold image and clean up segmentation with active contour

threshLevel = graythresh(fluorIm);
seg_im = im2bw(fluorIm,threshLevel);
seg_im = imfill(seg_im,'holes');
seg_im = medfilt2(seg_im,[5 5]);
sharpFluor = imsharpen(fluorIm,'radius',20,'amount',.1);
seg_im = activecontour(sharpFluor,imdilate(seg_im,strel('disk',10)), 50,'Chan-Vese');


%% collecting centroids via watershed transform on distance image
D = -bwdist(~seg_im);
D(~seg_im) = -Inf;
L = watershed(D);

% remove centroids with really small or really large areas
lc = regionprops(L,'Centroid','Area');
centroids = cat(1,lc.Centroid);
areas = cat(1,lc.Area);
centroids(areas>1e4 | areas<=50,:) = [];
centroids = round(centroids);

currentSeeds = centroids;

end
