
close all;
clear all;clc

CellLine = 'SUM149';
%load(fullfile(pwd, ['BootstrapFits_' CellLine '.mat']))
load(fullfile(pwd, ['BootstrapFits_IJcounts_' CellLine '.mat']))
%load(fullfile(pwd, ['BootstrapFits_IJcounts_' CellLine '_20161121.mat']))

c1 = [1 0 0];
plotLocs = [1 4 7 2 5 8 3 6 9];

LabelLocs_x = [400 400 -100 400 400 -100 400 -100 -100];
LabelLocs_y = [.2 .2 .9 .2 .2 .9 .2 .9 .9];
replicates = distinguishable_colors(7);
replicates(2,:) = [];

%%

wid = .1;
hei = .06;

allLocs = [.24 .717 wid hei;%col 1
    .24 .417 wid hei;
    .135 .257 wid hei;
    ...
    .52 .717 wid hei;%col 2
    .52 .417 wid hei;
    .415 .257 wid hei;
    ...
    .80 .717 wid hei;%col 3
    .695 .557 wid hei;
    .695 .257 wid hei;];



cf = figure(100);clf
set(cf,'Position',[1300 100 1400 800],...
    'Name', 'Figure 3: SUM149 Response')
cntr = 0;

for expTimes = 1:3
    for treatIter = [3 6 9]%2:length(uconcs)
        cntr = cntr+1;
        
        cFits = cat(2, bestFitModels_smooth{treatIter-1, :, expTimes});
        
        %ci = min(cellData.allYweights{expTimes}{treatIter});
        ci = floor(mean(cellData.allYweights{expTimes}{treatIter}));
        
        allTspan = cellData.Tspan{expTimes};
        postTreat_TP = cellData.DrugAdded_Tp + 1;
        pt_Tspan = allTspan(postTreat_TP:end);
        pt_Tspan_smooth = linspace(pt_Tspan(1), pt_Tspan(end), size(cFits,1));
        
        mtp = find(pt_Tspan_smooth >= allTspan(ci), 1);
        pt_Tspan_smooth = pt_Tspan_smooth(1:mtp);
        cFits = cFits(1:mtp,:);
        
        % confidence bound on model fits
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
        
        for wellIter = 1:6
            c_weights = ones(size(cellCounts));
            if isfield(cellData, 'allYweights')
                lvtp = cellData.allYweights{expTimes}{1}(wellIter);
                if lvtp < length(allTspan)
                    c_weights = zeros(size(cellCounts));
                    c_weights(1:(lvtp-cellData.DrugAdded_Tp)) = 1;
                end
            end
        end
        
        bnds = tinv([0.025, 0.975],nSamp-1)' * daStd/sqrt(nSamp);
        ub = bnds(2,:);
        lb = -bnds(1,:);
        
        subplot(3,3, plotLocs(cntr))
        errorbar(allTspan,daMean,...
            lb,ub,...
            'ks','MarkerSize',2,'LineWidth',2);
        hold on
%         % plot each replicate
%         for wi = 1:size(cellCounts,2)
%             cp = plot(allTspan(1:cellData.DrugAdded_Tp),...
%                 cellCounts(1:cellData.DrugAdded_Tp,wi),'k-');
%             set(cp,'Color',replicates(wi,:));
%             cp = plot(allTspan(cellData.DrugAdded_Tp+1:end),...
%                 cellCounts(cellData.DrugAdded_Tp+1:end,wi),'k-');
%             set(cp,'Color',replicates(wi,:));
%         end
        
        theAlphaVal = 0.1;
        fvad = .8.*ones(size(theShade,1),1);
        fvad(length(pt_Tspan_smooth)) = theAlphaVal;
        fvad(end) = theAlphaVal;
        
        %hold on
        %cp = patch(theShade(:,1),theShade(:,2),c1(1,:),'FaceAlpha',theAlphaVal,...
        %    'EdgeColor',c1(1,:),'EdgeAlpha','flat','LineStyle','-',...
        %    'FaceVertexAlphaData',fvad,'LineWidth',2,'LineStyle','-');
        %%set(currplot,'Color',c1(1,:))
        %hold off
        
        daTitle = sprintf('%3.0f nM',cellData.uconcs{expTimes}(treatIter));
        title(daTitle);
        
        ax = gca;
        set(gca,'XLim',[-120 (cellData.Tspan{expTimes}(end)+1*24)])
        set(gca,'YLim',[0 4e4])
        ax = gca;
        %set(gca,'YTickLabel',num2str(get(gca,'YTick').'))
        set(gca,'YTickLabel',sprintf('%2.1f\n',ax.YTick/1e4));
        fs = 12;
        text(-120, ax.YLim(2)*1.07, 'x10^4','FontSize',fs);
        
        cAUC = sprintf('%3.0f',cellData.uconcs{expTimes}(treatIter)*...
            cellData.ExposureTimes(expTimes));
        dox_pss = sprintf('%3.1f\n',cellData.dox_pss{expTimes}(treatIter));
        %text(LabelLocs_x(cntr), ax.YLim(2)*(LabelLocs_y(cntr)-.12), ['{\itAUC} = ' cAUC],'FontSize',fs)
        %text(LabelLocs_x(cntr), ax.YLim(2)*LabelLocs_y(cntr), ['{\itC_{B,max}} = ' dox_pss],'FontSize',fs)
        
        
        annLocs = allLocs(cntr,:);
        an = annotation('textbox', annLocs,'String',...
            ['{\itC_{B,max}} = ' dox_pss '{\itAUC} = ' cAUC],...
            'FontSize',fs,'FitBoxToText','off');
        set(an,'FontName','Times')
        
        get(gca,'FontSize');
        set(gca,'FontSize',fs);
    end
    
end
