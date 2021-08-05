function ys=mylowess(xy,xs,span)
%MYLOWESS Lowess smoothing, preserving x values
%   YS=MYLOWESS(XY,XS) returns the smoothed version of the x/y data in the
%   two-column matrix XY, but evaluates the smooth at XS and returns the
%   smoothed values in YS.  Any values outside the range of XY are taken to
%   be equal to the closest values.

if nargin<3 || isempty(span)
    span = .3;
end

if size(xy, 2) == 2
    % Sort and get smoothed version of xy data
    xy = sortrows(xy);
    x1 = xy(:,1);
    y1 = xy(:,2);
    ys1 = smooth(x1,y1,span,'loess');
    
    % Remove repeats so we can interpolate
    t = diff(x1)==0;
    x1(t)=[]; ys1(t) = [];
    
    % Interpolate to evaluate this at the xs values
    %ys = interp1(x1,ys1,xs,'pchip',NaN);
    ys = interp1(x1,ys1,xs);
    
    % Some of the original points may have x values outside the range of the
    % resampled data.  Those are now NaN because we could not interpolate them.
    % Replace NaN by the closest smoothed value.  This amounts to extending the
    % smooth curve using a horizontal line.
    if any(isnan(ys))
        ys(xs<x1(1)) = ys1(1);
        ys(xs>x1(end)) = ys1(end);
    end
    
elseif size(xy, 2) == 3
    
    % Normalize distances
    
    % Sort training values
    xy = sortrows(xy);
    x1 = xy(:,1:2);
    y1 = xy(:,end);
    
    % local regression to fit surface to points
    %sFit = fit(x1,y1,'lowess','Span',span);
    % evaluate fit at test points, xs
    %ys = sFit(xs);
    
    % use span to determine number of points to use for each test case
    numPts = ceil(size(x1,1) * span);
    if numPts < 3
        numPts = 3;
    end
    
    % initialize output vector
    ys = zeros(size(xs,1),1);
    
    % if points lie outsize training set, fix value to extremes
    ltp = xs(:,1) < x1(1,1) | xs(:,2) < x1(1,2);
    ys(ltp) = y1(1);
    gtp = xs(:,1) > x1(end,1) | xs(:,2) > x1(end,2);
    ys(gtp) = y1(end);
    
    % only caluculate lowess for points within training set bounds
    validPts = find(~(ltp | gtp));
    %validPts = find(true(size(ys)));
    
    % loop through valid points
    for vpi = 1:length(validPts)
        % calcluate distance to each point in test set
        pp = validPts(vpi);
        dadist = x1 - repmat(xs(pp,:), size(x1,1),1);
        dadist = sqrt(sum(dadist.^2, 2));
        
        [v,ioi] = sort(dadist);
        
        % calculate weighting vector
        maxD = v(numPts);
        w = (1 - (v(1:numPts)./maxD).^3).^3;
        
        % fit local linear regression (or moving average?!)
        ft = 'poly10';
        fo = fitoptions(ft);
        fo.Weights = w;
        sfit = fit(x1(ioi(1:numPts),:), y1(ioi(1:numPts)), ft,fo);
        
        
        ys(pp) = sfit(xs(pp,:));
        
        %testX = [linspace(min(x1(ioi(1:numPts),1)),max(x1(ioi(1:numPts),1)),100)',...
        %    linspace(min(x1(ioi(1:numPts),2)),max(x1(ioi(1:numPts),2)),100)'];
        %figure(100);clf;
        %subplot(1,2,1)
        %hold on
        %plot(x1(ioi(1:end),1), y1(ioi(1:end)),'k.')
        %plot(xs(pp,1), ys(pp),'bo')
        %plot(testX(:,1), sfit(testX),'k-')
        %subplot(1,2,2)
        %hold on
        %plot(x1(ioi(1:end),2), y1(ioi(1:end)),'k.')
        %plot(xs(pp,2), ys(pp),'bo')
        %plot(testX(:,2), sfit(testX),'k-')
        %pause
        
    end
    
    %figure(1);clf;
    %subplot(1,2,1)
    %hold on
    %plot(x1(:,1), y1,'k.')
    %plot(xs(:,1), ys,'bo')
    %subplot(1,2,2)
    %hold on
    %plot(x1(:,2), y1,'k.')
    %plot(xs(:,2), ys,'bo')
    
end
