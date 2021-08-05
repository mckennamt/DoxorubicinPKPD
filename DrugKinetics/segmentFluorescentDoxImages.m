
function [seg_im,currentSeeds,proceedToTL] = segmentFluorescentDoxImages(Im_cellFL,Im_aifFL, previousSeeds)

% get initial seeds from simple segmentation information
proceedToTL = false;

diffIm = Im_cellFL - Im_aifFL;
diffIm = diffIm + abs(min(diffIm(:)));
fluorIm = adapthisteq(diffIm./max(diffIm(:)),'NumTiles',[20 20],...
    'NBins',256,...
    'ClipLimit',.1,...
    'Distribution','rayleigh',...
    'Alpha',.5);
fluorIm = wiener2((fluorIm),[10 10]);

threshLevel = graythresh(fluorIm);
seg_im = im2bw(fluorIm,threshLevel);
seg_im = imfill(seg_im,'holes');
seg_im = medfilt2(seg_im,[5 5]);

sharpFluor = imsharpen(fluorIm,'radius',5,'amount',2);

seg_im = activecontour(sharpFluor,imdilate(seg_im,strel('disk',10)), 25,'edge');

% collecting centroids via watershed transform on distance image
D = -bwdist(~seg_im);
D(~seg_im) = -Inf;
L = watershed(D);

% remove centroids with really small or really large areas
lc = regionprops(L,'Centroid','Area');
centroids = cat(1,lc.Centroid);
areas = cat(1,lc.Area);
centroids(areas>1e4 | areas<=5,:) = [];
centroids = round(centroids);


%%view seeds
%fluorSeeds = zeros(size(Im_cellFL));
%seedInds = sub2ind(size(Im_cellFL),centroids(:,2),centroids(:,1));
%fluorSeeds(seedInds) = 1;
%fluorSeeds_dilated = imdilate(fluorSeeds,ones(30,30));
%lab = label2rgb(fluorSeeds_dilated,'hsv');
lab = label2rgb(seg_im,'hsv');
figure(1);clf
imagesc(diffIm)
%caxis([0 20])
colormap gray
hold on
h = imshow(lab);
set(h,'AlphaData',.1)

istr = sprintf(['Good Segmentation? Press return if segmentation is good.',...
    '\nPress 0 if want to modify seeds.\nPress 1 to clear seeds and use TL.',...
    '\nPress 2 to use TL seeds on fluorescent images: ']);
goodSegmentation = input(istr);

if ~isempty(goodSegmentation)
    if goodSegmentation == 0
        [xb,yb] = modifySegmentationSeeds(fluorIm, centroids);
        seg_im = manualSeeds_fluorIm(fluorIm,xb, yb);
        centroids = round([xb,yb]);
        currentSeeds = centroids;
    elseif goodSegmentation == 1
        currentSeeds = previousSeeds;
        proceedToTL = true;
        return
    elseif goodSegmentation == 2
        [xb,yb] = modifySegmentationSeeds(fluorIm, previousSeeds);
        seg_im = manualSeeds_fluorIm(fluorIm,xb, yb);
        centroids = round([xb,yb]);
        currentSeeds = centroids;
    end
end

currentSeeds = centroids;

% use gradients from intensity image to refine segmentation

% fluorIm = medfilt2(fluorIm,[10 10]);
% diff_im = vuCurvatureAnisotropicDiffusion(fluorIm,'numIterations',10);%100
% grad_im = vuGradientMagnitude(diff_im);
%
% sig_im = vuSigmoidFunction(grad_im,'alpha',5,'beta',50);%10, 20
%
% %figure(1);clf;subplot(121);imagesc(sig_im);caxis([0 2])
% %subplot(122);imagesc(fluorIm)
%
% seg_im = vuFastMarchingLevelSet(sig_im,currentSeeds,'stoptime',20);%200
%
% threshLevel = graythresh(seg_im./max(seg_im(:)));
% seg_im = imcomplement(im2bw(seg_im./max(seg_im(:)),threshLevel));
% seg_im = imfill(seg_im,'holes');

%seg_im = activecontour(fluorIm,seg_im, 25,'edge');
%bw2 = activecontour(fluorIm,initMask, 50,'edge');

%figure(12);clf
%lab = label2rgb(seg_im,'autumn');
%imagesc(fluorIm)
%colormap gray
%hold on
%h = imshow(lab);
%set(h,'AlphaData',.1)
end

function [seg_im] = manualSeeds_fluorIm(fluorIm,xb, yb)

    seedInds = sub2ind(size(fluorIm),floor(yb),floor(xb));
    seedMask = zeros(size(fluorIm));
    seedMask(seedInds) = 1;
    dilStrel = strel('disk',25);
    seedMask = imdilate(seedMask,dilStrel);
    seg_im = activecontour(fluorIm,seedMask, 50,'edge');
    
end
