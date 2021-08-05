



CellLine = 'SUM149';
load(fullfile(pwd, [CellLine '_Predictions']))

cf = figure(121);clf
c1 = [1 0 0];

cntr = 0;
expTimes = 2;%index to select [6 12 24] hrs


for treatIter = [5 7 9]%2:length(uconcs)
    cntr = cntr+1;
    
    cFits = cat(2, predictedTimecourses_smooth{treatIter-1, :, expTimes});
    
    allTspan = cellData.Tspan{expTimes};
    postTreat_TP = cellData.DrugAdded_Tp + 1;
    pt_Tspan = allTspan(postTreat_TP:end);
    pt_Tspan_smooth = linspace(pt_Tspan(1), pt_Tspan(end), size(cFits,1));
    
    inds_95CI = ceil(.025*size(cFits,2));
    theShade = zeros(length(pt_Tspan_smooth)*2,2);
    theShade(:,1) = [pt_Tspan_smooth'; pt_Tspan_smooth(end:-1:1)'];
    for makeShade = 1:length(pt_Tspan_smooth)
        srted = sort(cFits(makeShade,:));
        theShade(makeShade,2) = srted(inds_95CI);
    end
    
    inds_95CI = floor(.975*size(cFits,2));
    shadeIndex = length(pt_Tspan_smooth);
    for makeShade = length(pt_Tspan_smooth):-1:1
        shadeIndex = shadeIndex+1;
        srted = sort(cFits(makeShade,:));
        theShade(shadeIndex,2) = srted(inds_95CI);
    end
    
    
    cellCounts = cellData.allYall{expTimes}{treatIter};
    daStd = std(cellCounts');
    daMean = mean(cellCounts');
    nSamp = size(cellCounts,2);
    
    bnds = tinv([0.025, 0.975],nSamp-1)' * daStd/sqrt(nSamp);
    ub = bnds(2,:);
    lb = -bnds(1,:);
    
    subplot(1,3, cntr)
    errorbar(allTspan,daMean,...
        lb,ub,...
        'ks','MarkerSize',2,'LineWidth',2)
    
    fvad = .8.*ones(size(theShade,1),1);
    fvad(length(pt_Tspan_smooth)) = .2;
    fvad(end) = .2;
    
    hold on
    theAlphaVal = .1;
    patch(theShade(:,1),theShade(:,2),c1(1,:),'FaceAlpha',theAlphaVal,...
        'EdgeColor',c1(1,:),'EdgeAlpha','flat','LineStyle','-',...
        'FaceVertexAlphaData',fvad)
    %set(currplot,'Color',c1(1,:))
    hold off
    
    daTitle = sprintf('%3.0f nM',cellData.uconcs{expTimes}(treatIter));
    title(daTitle)
    
    ax = gca;
    set(gca,'XLim',[-120 (allTspan(end)+1*24)])
    set(gca,'YLim',[0 2.5e4])
    %set(gca,'YTickLabel',num2str(get(gca,'YTick').'))
    set(gca,'YTickLabel',sprintf('%2.1f\n',ax.YTick/1e4))
    text(-120, 2.5e4*1.04, 'x10^4','FontSize',ax.FontSize)
end

set(cf,'Position',[1100 200 1600 350])

%% make single prediction, top right panel

CellLine = 'SUM149';
load(fullfile(pwd, [CellLine '_Predictions']))

cf = figure(122);clf
c1 = [1 0 0];

expTimes = 2;%index to select [6 12 24] hrs


for treatIter = 6%2:length(uconcs)
    
    cFits = cat(2, predictedTimecourses_smooth{treatIter-1, :, expTimes});
    
    allTspan = cellData.Tspan{expTimes};
    postTreat_TP = cellData.DrugAdded_Tp + 1;
    pt_Tspan = allTspan(postTreat_TP:end);
    pt_Tspan_smooth = linspace(pt_Tspan(1), pt_Tspan(end), size(cFits,1));
    
    inds_95CI = ceil(.025*size(cFits,2));
    theShade = zeros(length(pt_Tspan_smooth)*2,2);
    theShade(:,1) = [pt_Tspan_smooth'; pt_Tspan_smooth(end:-1:1)'];
    for makeShade = 1:length(pt_Tspan_smooth)
        srted = sort(cFits(makeShade,:));
        theShade(makeShade,2) = srted(inds_95CI);
    end
    
    inds_95CI = floor(.975*size(cFits,2));
    shadeIndex = length(pt_Tspan_smooth);
    for makeShade = length(pt_Tspan_smooth):-1:1
        shadeIndex = shadeIndex+1;
        srted = sort(cFits(makeShade,:));
        theShade(shadeIndex,2) = srted(inds_95CI);
    end
    
    
    cellCounts = cellData.allYall{expTimes}{treatIter};
    daStd = std(cellCounts');
    daMean = mean(cellCounts');
    nSamp = size(cellCounts,2);
    
    bnds = tinv([0.025, 0.975],nSamp-1)' * daStd/sqrt(nSamp);
    ub = bnds(2,:);
    lb = -bnds(1,:);
    
    errorbar(allTspan,daMean,...
        lb,ub,...
        'ks','MarkerSize',2,'LineWidth',2)
    
    fvad = .8.*ones(size(theShade,1),1);
    fvad(length(pt_Tspan_smooth)) = .2;
    fvad(end) = .2;
    
    hold on
    theAlphaVal = .1;
    patch(theShade(:,1),theShade(:,2),c1(1,:),'FaceAlpha',theAlphaVal,...
        'EdgeColor',c1(1,:),'EdgeAlpha','flat','LineStyle','-',...
        'FaceVertexAlphaData',fvad)
    %set(currplot,'Color',c1(1,:))
    hold off
    
    daTitle = sprintf('%3.0f nM',cellData.uconcs{expTimes}(treatIter));
    title(daTitle)
    
    ax = gca;
    set(gca,'XLim',[-120 (allTspan(end)+1*24)])
    set(gca,'YLim',[0 2.5e4])
    set(gca,'FontSize',14)
    %set(gca,'YTickLabel',num2str(get(gca,'YTick').'))
    set(gca,'YTickLabel',sprintf('%2.1f\n',ax.YTick/1e4))
    text(-120, 2.5e4*1.04, 'x10^4','FontSize',ax.FontSize)
end

set(cf,'Position',[100 100 800 500])