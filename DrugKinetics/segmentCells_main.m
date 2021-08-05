

function [seg_im,currentSeeds, croppedCell] = ...
    segmentCells_main(Im_cellTL, Im_cellFL, Im_aifFL, ...
    probGoodFluorescentSignal, numCells, previousSeeds, Im_cellTL_old,...
    headlessSegmentation, headlessSeeds, croppedCell)

%% Segment cells for doxorubicin uptake studies


if probGoodFluorescentSignal
    
    %[seg_im,currentSeeds,proceedToTL] = segmentFluorescentDoxImages(Im_cellFL,Im_aifFL,previousSeeds);
    [seg_im,currentSeeds,proceedToTL, reviewSegmentation] = ...
        segmentFluorescentDoxImages_boosted(Im_cellFL,Im_aifFL,previousSeeds,...
        headlessSegmentation, headlessSeeds, Im_cellTL);
    
    if proceedToTL
        needToSegment = true;
    elseif reviewSegmentation
        figure(2);clf
        lab = label2rgb(seg_im,'autumn');
        subplot(1,2,1)
        imagesc(Im_cellTL)
        colormap gray
        hold on
        h = imshow(lab);
        set(h,'AlphaData',.2)
        
        subplot(1,2,2)
        imagesc(Im_cellFL)
        colormap gray
        hold on
        h = imshow(lab);
        set(h,'AlphaData',.2)
        
        istr = sprintf('Good Segmentation? Press return if segmentation is good.\nPress 0 if want to modify seeds.\nPress 1 to clear seeds: ');
        goodSegmentation = input(istr);
        
        if isempty(goodSegmentation)
            needToSegment = false;
        elseif goodSegmentation == 0
            needToSegment = true;
            previousSeeds = currentSeeds;
        elseif goodSegmentation == 1
            needToSegment = true;
        end
    else
        needToSegment = false;
    end
    
else
    needToSegment = true;
end

while needToSegment
    %% get seeds for segmentation
    if true%isempty(previousSeeds)
        
        [previousSeeds,croppedCell] = quickSeedGeneration(Im_cellTL,numCells,croppedCell);
        
        %istr = sprintf('Update seeds? Press 1, otherwise enter: ');
        %goodSegmentation = input(istr);
        %
        %if isempty(goodSegmentation) && ~isempty(previousSeeds)
        %    needToSegment = false;
        %elseif ~isempty(goodSegmentation) || isempty(previousSeeds)
        %    [previousSeeds,croppedCell] = quickSeedGeneration(Im_cellTL,numCells,croppedCell);
        %end
    else
        % in case image shifted, use x-corr to find shift
        %previousSeeds = seedShift(Im_cellTL_old, Im_cellTL, previousSeeds);
    end
    
    
    % manually add in or remove seeds
    [xb,yb] = modifySegmentationSeeds(Im_cellTL, previousSeeds);
    
    currentSeeds = [xb yb];
    
    %% segement image
    
    % use seeds to segment image
    %diff_im = vuCurvatureAnisotropicDiffusion(Im_cellTL,'numIterations',10);
    %grad_im = vuGradientMagnitude(diff_im);
    %
    %sig_im = vuSigmoidFunction(grad_im,'alpha',5,'beta',20);
    %%figure(121);imagesc(sig_im)
    %seg_im = vuFastMarchingLevelSet(sig_im,[xb yb],'stoptime',100);
    %
    %threshLevel = graythresh(seg_im./max(seg_im(:)));
    %seg_im = imcomplement(im2bw(seg_im./max(seg_im(:)),threshLevel));
    %seg_im = imfill(seg_im,'holes');
    
    seedInds = sub2ind(size(Im_cellTL),floor(yb),floor(xb));
    seedMask = zeros(size(Im_cellTL));
    seedMask(seedInds) = 1;
    dilStrel = strel('disk',25);
    seedMask = imdilate(seedMask,dilStrel);
    
    
    modIm = adapthisteq(Im_cellTL./max(Im_cellTL(:)),'NumTiles',[20 20],...
        'Distribution','rayleigh','Alpha',.01);
    %modIm = histeq(modIm,12);
    modIm = wiener2(modIm,[10 10]);
    modIm = imadjust(modIm);
    %keyboard
    seg_im = activecontour(modIm,seedMask, 50,'edge');
    
    
    %figure(1);clf;
    %imagesc(Im_cellTL);colormap gray
    %title('Segmented Cells')
    %hold on
    %lab = label2rgb(seg_im,'autumn');
    %h = imshow(lab);
    %set(h,'AlphaData',.2);
    
    figure(2);clf
    lab = label2rgb(seg_im,'autumn');
    subplot(1,2,1)
    imagesc(Im_cellTL)
    colormap gray
    hold on
    h = imshow(lab);
    set(h,'AlphaData',.2)
    
    subplot(1,2,2)
    imagesc(Im_cellFL)
    colormap gray
    hold on
    h = imshow(lab);
    set(h,'AlphaData',.2)
    
    
    goodSegmentation = input('Good Segmentation? Press ''return'' if segmentation is good. Otherwise press 0: ');
    
    if isempty(goodSegmentation)
        needToSegment = false;
        
        % collecting centroids via watershed transform on distance image
        D = -bwdist(~seg_im);
        D(~seg_im) = -Inf;
        L = watershed(D);
        
        % remove centroids with really small or really large areas
        lc = regionprops(L,'Centroid','Area');
        centroids = cat(1,lc.Centroid);
        areas = cat(1,lc.Area);
        centroids(areas>1e4 | areas<=500,:) = [];
        %currentSeeds = round(centroids);
        
    else
        needToSegment = true;
        previousSeeds = currentSeeds;
    end
end

% D = bwdist(~seg_im);
% D = -D;
% D(~seg_im) = -Inf;
% L = watershed(D);
%
% % 'refresh' centroid list to reflect true position of cells
% lc = regionprops(L,'Centroid');
% centroids = cat(1,lc.Centroid);
% currentSeeds = round(centroids);
