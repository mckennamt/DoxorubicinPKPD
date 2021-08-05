

function [xb,yb] = modifySegmentationSeeds(Im2Segment, previousSeeds)

prevSeeds = zeros(size(Im2Segment));
seedInds = sub2ind(size(Im2Segment),floor(previousSeeds(:,2)),floor(previousSeeds(:,1)));
prevSeeds(seedInds) = 1;
prevSeeds_dilated = imdilate(prevSeeds,ones(30,30));

xb = previousSeeds(:,1);
yb = previousSeeds(:,2);

figure(1);clf;
imagesc(Im2Segment)
%caxis([650 1880])
%caxis([1800 3500])
hold on
h = imshow(cat(3,ones(size(Im2Segment)),zeros(size(Im2Segment)),zeros(size(Im2Segment))));
set(h,'AlphaData',.2.*prevSeeds_dilated);

colormap gray
title('Select seeds for cells')
cbutton = 1;


while ~isempty(cbutton)
    
    set(h,'AlphaData',.2.*prevSeeds_dilated);
    
    [x,y, cbutton] = ginput(1);
    
    if isempty(xb) && ~isempty(x)
        xb = [xb;x];
        yb = [yb;y];
        seedInds = sub2ind(size(Im2Segment),floor(y),floor(x));
        prevSeeds(seedInds) = 1;
        prevSeeds_dilated(floor(y)-15:floor(y)+15,floor(x)-15:floor(x)+15) = 1;
        continue
    end
    
    if ~isempty(x)
        distToPrev = sqrt(sum(([xb yb] - repmat([x,y], length(xb),1)).^2,2));
        removeMe = find(distToPrev < 15,1);
        
        if sum(removeMe)>0
            
            seedInds = sub2ind(size(Im2Segment),floor(yb(removeMe)),floor(xb(removeMe)));
            prevSeeds(seedInds) = 1;
            
            y_min = max([floor(yb(removeMe))-15 1]);
            y_max = min([floor(yb(removeMe))+15 size(prevSeeds_dilated,1)]);
            
            x_min = max([floor(xb(removeMe))-15 1]);
            x_max = min([floor(xb(removeMe))+15 size(prevSeeds_dilated,2)]);
            
            prevSeeds_dilated(y_min:y_max,x_min:x_max) = 0;
            
            xb(removeMe) = [];
            yb(removeMe) = [];
            
        else
            
            xb = [xb;x];
            yb = [yb;y];
            
            seedInds = sub2ind(size(Im2Segment),floor(y),floor(x));
            prevSeeds(seedInds) = 1;
            
            y_min = max([floor(y)-15 1]);
            y_max = min([floor(y)+15 size(prevSeeds_dilated,1)]);
            
            x_min = max([floor(x)-15 1]);
            x_max = min([floor(x)+15 size(prevSeeds_dilated,2)]);
            
            prevSeeds_dilated(y_min:y_max,x_min:x_max) = 1;
            
        end
    end
end