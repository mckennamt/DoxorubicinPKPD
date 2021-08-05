
close all;
clear all;clc

cf = figure(102);clf
set(cf,'Position',[1300 100 1400 800],...
    'Name','Figure 4: All Cell Line Response')

allCellLine = {'SUM149','MDAMB231','MDAMB453','MDAMB468'};

plotLocs = [1 5 9 2 6 10 3 7 11 4 8 12];

c1 = [1 0 0];

expTimes = 1;%only 6 hour dataset
drugSel = [6 7 8];
%drugSel = [2 3 4];

yAxisLims = [3.4e4 1.8e4 2.5e4 3e4];
LabelLocs_x = [300, 300, -100,...
    -100, -100, -100,...
    -100, -100, -100,...
    -100, -100, -100];
LabelLocs_y = [.2, .2, .9,...
    .9, .9, .9,...
    .9, .9, .9,...
    .9, .9, .9];
replicates = distinguishable_colors(7);
replicates(2,:) = [];


tr = .858;
mr = .557;
br = .257;

wid = .1;%089
hei = .06;

allLocs = [.185 .717 wid hei;%col 1
    .185 .417 wid hei;
    .135 br wid hei;
    ...
    .34 tr wid hei;%col 2
    .34 mr wid hei;
    .34 br wid hei;
    ...
    .545 tr wid hei;%col 3
    .545 mr wid hei;
    .545 br wid hei;
    ...
    .753 tr wid hei;%col 4
    .753 mr wid hei;
    .753 br wid hei;];


cntr = 0;

for cellLineSelect = 1:4
    
    CellLine = allCellLine{cellLineSelect};
    %load(fullfile(pwd, ['BootstrapFits_' CellLine '.mat']))
    load(fullfile(pwd, ['BootstrapFits_IJcounts_' CellLine '.mat']))
    %load(fullfile(pwd, ['BootstrapFits_IJcounts_' CellLine '_20161121.mat']))
    
    kp
    
    for treatIter = [6 7 8]%[6 7 8]
        cntr = cntr+1;
        
        cellData.allYweights{expTimes}{treatIter};
        %ci = min(cellData.allYweights{expTimes}{treatIter});
        ci = floor(mean(cellData.allYweights{expTimes}{treatIter}));
        
        cFits = cat(2, bestFitModels_smooth{treatIter-1, :, expTimes});
        allTspan = cellData.Tspan{expTimes};
        postTreat_TP = cellData.DrugAdded_Tp + 1;
        pt_Tspan = allTspan(postTreat_TP:end);
        pt_Tspan_smooth = linspace(pt_Tspan(1), pt_Tspan(end), size(cFits,1));
        
        mtp = find(pt_Tspan_smooth >= allTspan(ci), 1);
        if mtp < length(pt_Tspan_smooth)
            mtp = mtp+1;
        end
        pt_Tspan_smooth = pt_Tspan_smooth(1:mtp);
        cFits = cFits(1:mtp,:);
        
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
        
        % show cell counts
        %ci = min(cellData.allYweights{expTimes}{treatIter});
        cellCounts = cellData.allYall{expTimes}{treatIter}(1:ci,:);
        allTspan = allTspan(1:ci);
        daStd = std(cellCounts');
        daMean = mean(cellCounts');
        nSamp = size(cellCounts,2);
        
        bnds = tinv([0.025, 0.975],nSamp-1)' * daStd/sqrt(nSamp);
        ub = bnds(2,:);
        lb = -bnds(1,:);
        
        subplot(3,4, plotLocs(cntr))
        errorbar(allTspan,daMean,...
            lb,ub,...
            'ks','MarkerSize',2,'LineWidth',2)
        hold on
        % plot each replicate
%         for wi = 1:size(cellCounts,2)
%            cp = plot(allTspan(1:cellData.DrugAdded_Tp),...
%                cellCounts(1:cellData.DrugAdded_Tp,wi),'k-');
%            set(cp,'Color',replicates(wi,:));
%            cp = plot(allTspan(cellData.DrugAdded_Tp+1:end),...
%                cellCounts(cellData.DrugAdded_Tp+1:end,wi),'k-');
%            set(cp,'Color',replicates(wi,:));
%         end
        
        theAlphaVal = 0.1;
        fvad = .8.*ones(size(theShade,1),1);
        fvad(length(pt_Tspan_smooth)) = theAlphaVal;
        fvad(end) = theAlphaVal;
        
        hold on
        patch(theShade(:,1),theShade(:,2),c1(1,:),'FaceAlpha',theAlphaVal,...
            'EdgeColor',c1(1,:),'EdgeAlpha','flat','LineStyle','-',...
            'FaceVertexAlphaData',fvad,'LineWidth',2,'LineStyle','-');
        %set(currplot,'Color',c1(1,:))
        hold off
        
        daTitle = sprintf('%3.0f nM',cellData.uconcs{expTimes}(treatIter));
        title(daTitle)
        
        set(gca,'XLim',[-120 (cellData.Tspan{expTimes}(end)+1*24)])
        %set(gca,'XLim',[-120 600])
        set(gca,'YLim',[0 yAxisLims(cellLineSelect)])
        
        %set(gca,'YTickLabel',num2str(get(gca,'YTick').'))
        ax = gca;
        %set(gca,'YTick',ax.YTick)
        set(gca,'YTickLabel',sprintf('%2.1f\n',ax.YTick/1e4))
        fs = 12;
        text(-120, yAxisLims(cellLineSelect)*1.07, 'x10^4','FontSize',fs)
        set(gca,'FontSize',fs)
        
        cAUC = sprintf('%3.0f',cellData.uconcs{expTimes}(treatIter)*...
            cellData.ExposureTimes(expTimes));
        dox_pss = sprintf('%3.1f\n',cellData.dox_pss{expTimes}(treatIter));
        %text(LabelLocs_x(cntr), ax.YLim(2)*(LabelLocs_y(cntr)-.12), ['{\itAUC} = ' cAUC],'FontSize',fs)
        %text(LabelLocs_x(cntr), ax.YLim(2)*LabelLocs_y(cntr), ['{\itC_{B,max}} = ' dox_pss],'FontSize',fs)
        %set(gca,'YTickLabelMode','Auto')
        
        annLocs = allLocs(cntr,:);
        an = annotation('textbox', annLocs,'String',...
            ['{\itC_{B,max}} = ' dox_pss '{\itAUC} = ' cAUC],...
            'FontSize',fs,'FitBoxToText','off','BackgroundColor','w',...
            'FaceAlpha',.8);
        set(an,'FontName','Times')
        
    end
    
end
