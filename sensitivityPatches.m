

function figHandle = sensitivityPatches(figHandle, currAx,...
    holdSens, l_cb_smoove,...
    patchFaceAlpha, patchEdgeAlpha, sensThresh)

% add gray region to denote parameter insensitive over range

%sf = currAx.YLim(2) - currAx.YLim(1);
%semilogx(l_cb_smoove, (smooth(holdSens)*sf + currAx.YLim(1)),'b--');

ioi = find(smooth(holdSens) < sensThresh);

if isempty(ioi)
    return
end

segments_co = find(diff(ioi)>1);
segments = cell(length(segments_co) + 1, 1);

if isempty(segments_co)
    segments{1} = ioi;
    
    
else
    
    segments{1} = ioi(1:segments_co(1));
    segIter = 1;
    while segIter <= length(segments_co)
        if segIter+1 > length(segments_co)
            segments{segIter+1} = ioi(segments_co(segIter)+1:end);
        else
            segments{segIter+1} = ioi(segments_co(segIter)+1:segments_co(segIter+1));
        end
        segIter = segIter + 1;
    end
end


for iter = 1:length(segments)
    
    c_segment = segments{iter};
    
    x_sens_1 = l_cb_smoove(c_segment(1));
    x_sens_2 = l_cb_smoove(c_segment(end));
    
    if c_segment(1) == 1
        x_sens_1 = currAx.XLim(1);
    elseif c_segment(end) == length(holdSens)
        x_sens_2 = currAx.XLim(2);
    end
    
    semilogx(ones(1,100).*x_sens_1, linspace(currAx.YLim(1), currAx.YLim(2)),'k-')
    semilogx(ones(1,100).*x_sens_2, linspace(currAx.YLim(1), currAx.YLim(2)),'k-')
    pv = patch([x_sens_1 x_sens_2 x_sens_2 x_sens_1]',...
        [currAx.YLim(1) currAx.YLim(1) currAx.YLim(2) currAx.YLim(2)]','k');
    set(pv,'FaceAlpha',patchFaceAlpha,'EdgeAlpha',patchEdgeAlpha)
    
end