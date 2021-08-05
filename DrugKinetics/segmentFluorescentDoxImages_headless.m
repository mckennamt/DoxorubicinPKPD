
function [seg_im,currentSeeds,proceedToTL] = segmentFluorescentDoxImages_headless(Im_cellFL,Im_aifFL)

% get initial seeds from simple segmentation information
proceedToTL = false;

% subtract off AIF function
diffIm = Im_cellFL - Im_aifFL;
diffIm = diffIm + abs(min(diffIm(:)));

% normalize image intensities
fluorIm = adapthisteq(diffIm./max(diffIm(:)),'NumTiles',[20 20],...
    'NBins',256,...
    'ClipLimit',.1,...
    'Distribution','rayleigh',...
    'Alpha',.5);

% subtract off background with rolling ball filter
diskSize = 20;
rbb = imfilter(fluorIm, fspecial('disk', diskSize), 'replicate','same');
fluorIm = fluorIm - rbb;

% clean up noise
fluorIm = wiener2((fluorIm),[5 5]);

%% threshold image
% for all other cell lines
threshLevel = graythresh(fluorIm);
%%for MDAMB231
%threshLevel2 = multithresh(fluorIm,3);

seg_im = im2bw(fluorIm,threshLevel);
%seg_im = im2bw(fluorIm,threshLevel2(3));
%seg_im = imfill(seg_im,'holes');

seg_im = medfilt2(seg_im,[10 10]);

% remove small areas
seg_im = bwareaopen(seg_im,100);

%sharpFluor = imsharpen(fluorIm,'radius',5,'amount',5);
sharpFluor = fluorIm;

seg_im = activecontour(sharpFluor,imdilate(seg_im,strel('disk',7)), 10,'edge');

%%
% collecting centroids via watershed transform on distance image
D = -bwdist(~seg_im);
D(~seg_im) = -Inf;
L = watershed(D);

% remove centroids with really small or really large areas
lc = regionprops(L,'Centroid','Area');
centroids = cat(1,lc.Centroid);
areas = cat(1,lc.Area);
centroids(areas>1e4 | areas<=100,:) = [];
centroids = round(centroids);


currentSeeds = centroids;


end
