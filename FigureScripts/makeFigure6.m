
close all
clear all;clc

CellLine = 'SUM149';
%load(fullfile(pwd, [CellLine '_Predictions']))
load([CellLine '_Predictions_IJCounts_only12.mat']);
% %%
% cf = figure(6);clf
% set(cf,'Position',[1450 250 1250 500],...
%     'Name', 'Figure 6: Response Prediction')
% fs = 12;
% c1 = [1 0 0];
% 
% % Concentrations = [0 10 20 39 78 156 313 625 1250 2500]
% treatIters = [6 7 9 4 6 7];
% expTimes = [1 1 1 2 2 2];%index to select [6 24] hrs test groups
% 
% wid = .11;
% hei = .1;
% allLocs = [.23 .6 wid hei;%col 1
%     .23 .12 wid hei;
%     ...
%     .42 .81 wid hei;%col 2
%     .42 .33 wid hei;
%     ...
%     .70 .81 wid hei;%col 3
%     .70 .33 wid hei];
% allLoc_loc = [1 3 5 2 4 6];
% 
% cntr = 0;
% %expTime = 2;%index to select [6 12 24] hrs
% 
% for spi = 1:6
%     
%     treatIter = treatIters(spi);
%     expTime = expTimes(spi);
%     
%     cFits = cat(2, predictedTimecourses_smooth{treatIter-1, :, expTime});
%     
%     allTspan = cellData.Tspan{expTime};
%     postTreat_TP = cellData.DrugAdded_Tp + 1;
%     pt_Tspan = allTspan(postTreat_TP:end);
%     pt_Tspan_smooth = linspace(pt_Tspan(1), pt_Tspan(end), size(cFits,1));
%     
%     inds_95CI = ceil(.025*size(cFits,2));
%     theShade = zeros(length(pt_Tspan_smooth)*2,2);
%     theShade(:,1) = [pt_Tspan_smooth'; pt_Tspan_smooth(end:-1:1)'];
%     for makeShade = 1:length(pt_Tspan_smooth)
%         srted = sort(cFits(makeShade,:));
%         theShade(makeShade,2) = srted(inds_95CI);
%     end
%     
%     inds_95CI = floor(.975*size(cFits,2));
%     shadeIndex = length(pt_Tspan_smooth);
%     for makeShade = length(pt_Tspan_smooth):-1:1
%         shadeIndex = shadeIndex+1;
%         srted = sort(cFits(makeShade,:));
%         theShade(shadeIndex,2) = srted(inds_95CI);
%     end
%     
%     
%     % show cell counts
%     ci = min(cellData.allYweights{expTime}{treatIter});
%     cellCounts = cellData.allYall{expTime}{treatIter}(1:ci,:);
%     allTspan = allTspan(1:ci);
%     daStd = std(cellCounts');
%     daMean = mean(cellCounts');
%     nSamp = size(cellCounts,2);
%     
%     for wellIter = 1:6
%         c_weights = ones(size(cellCounts));
%         if isfield(cellData, 'allYweights')
%             lvtp = cellData.allYweights{expTime}{1}(wellIter);
%             if lvtp < length(allTspan)
%                 c_weights = zeros(size(cellCounts));
%                 c_weights(1:(lvtp-cellData.DrugAdded_Tp)) = 1;
%             end
%         end
%     end
%     
%     bnds = tinv([0.025, 0.975],nSamp-1)' * daStd/sqrt(nSamp);
%     ub = bnds(2,:);
%     lb = -bnds(1,:);
%     
% %     cellCounts = cellData.allYall{expTime}{treatIter};
% %     daStd = std(cellCounts');
% %     daMean = mean(cellCounts');
% %     nSamp = size(cellCounts,2);
% %     
% %     bnds = tinv([0.025, 0.975],nSamp-1)' * daStd/sqrt(nSamp);
% %     ub = bnds(2,:);
% %     lb = -bnds(1,:);
%     
%     subplot(2,3, spi)
%     errorbar(allTspan,daMean,...
%         lb,ub,...
%         'ks','MarkerSize',2,'LineWidth',2)
%     
%     theAlphaVal = 0.1;
%     fvad = .8.*ones(size(theShade,1),1);
%     fvad(length(pt_Tspan_smooth)) = theAlphaVal;
%     fvad(end) = theAlphaVal;
%     
%     hold on
%     patch(theShade(:,1),theShade(:,2),c1(1,:),'FaceAlpha',theAlphaVal,...
%         'EdgeColor',c1(1,:),'EdgeAlpha','flat','LineStyle','-',...
%         'FaceVertexAlphaData',fvad)
%     %set(currplot,'Color',c1(1,:))
%     hold off
%     
%     daTitle = sprintf('%3.0f nM',cellData.uconcs{expTime}(treatIter));
%     title(daTitle)
%     
%     cAUC = sprintf('%3.0f',cellData.uconcs{expTime}(treatIter)*...
%         cellData.ExposureTimes(expTime));
%     dox_pss = sprintf('%3.1f\n',cellData.dox_pss{expTime}(treatIter));
%     %text(LabelLocs_x(cntr), ax.YLim(2)*(LabelLocs_y(cntr)-.12), ['{\itAUC} = ' cAUC],'FontSize',fs)
%     %text(LabelLocs_x(cntr), ax.YLim(2)*LabelLocs_y(cntr), ['{\itC_{B,max}} = ' dox_pss],'FontSize',fs)
%     
%     
%     annLocs = allLocs(allLoc_loc(spi),:);
%     an = annotation('textbox', annLocs,'String',...
%         ['{\itC_{B,max}} = ' dox_pss '{\itAUC} = ' cAUC],...
%         'FontSize',fs,'FitBoxToText','off');
%     set(an,'FontName','Times')
%     
%     ax = gca;
%     set(gca,'XLim',[-100 (cellData.Tspan{expTime}(end)+1*24)])
%     set(gca,'YLim',[0 4e4])
%     
%     set(gca,'YTickLabel',sprintf('%2.1f\n',ax.YTick/1e4))
%     set(gca,'FontSize',fs)
%     text(-100, ax.YLim(2)*1.07, 'x10^4','FontSize',fs)
% end
% 

%%

cf = figure(66);clf
set(cf,'Position',[1300 100 800 800],...
    'Name', 'Figure 6: Response Prediction')
fs = 12;
c1 = [1 0 0];

treatIters = [7 8 9 6 7 10];
expTimes = [1 1 1 3 3 3];%index to select [6 12 24] hrs
pred_expTimes = [1 1 1 2 2 2];%index to select [6 12 24] hrs

plotLocs = [1 3 5 2 4 6];

wid = .17;
hei = .06;
allLocs = [.29 .72 wid hei;%col 1
    .14 .55 wid hei;
    .14 .25 wid hei;
    ...
    .58 .86 wid hei;%col 2
    .58 .55 wid hei;
    .58 .25 wid hei];

allLoc_loc = [1 2 3 4 5 6];

cntr = 0;
%expTime = 2;%index to select [6 12 24] hrs

for spi = 1:6

    treatIter = treatIters(spi);
    expTime = expTimes(spi);
    p_expTime = pred_expTimes(spi);

    cFits = cat(2, predictedTimecourses_smooth{treatIter-1, :, p_expTime});

    allTspan = cellData.Tspan{expTime};
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


    % show cell counts
    ci = min(cellData.allYweights{expTime}{treatIter});
    cellCounts = cellData.allYall{expTime}{treatIter}(1:ci,:);
    allTspan = allTspan(1:ci);
    daStd = std(cellCounts');
    daMean = mean(cellCounts');
    nSamp = size(cellCounts,2);
    
    for wellIter = 1:6
        c_weights = ones(size(cellCounts));
        if isfield(cellData, 'allYweights')
            lvtp = cellData.allYweights{expTime}{1}(wellIter);
            if lvtp < length(allTspan)
                c_weights = zeros(size(cellCounts));
                c_weights(1:(lvtp-cellData.DrugAdded_Tp)) = 1;
            end
        end
    end
    
    bnds = tinv([0.025, 0.975],nSamp-1)' * daStd/sqrt(nSamp);
    ub = bnds(2,:);
    lb = -bnds(1,:);
    
%     cellCounts = cellData.allYall{expTime}{treatIter};
%     daStd = std(cellCounts');
%     daMean = mean(cellCounts');
%     nSamp = size(cellCounts,2);
%     
%     bnds = tinv([0.025, 0.975],nSamp-1)' * daStd/sqrt(nSamp);
%     ub = bnds(2,:);
%     lb = -bnds(1,:);

    subplot(3,2, plotLocs(spi))
    errorbar(allTspan,daMean,...
        lb,ub,...
        'ks','MarkerSize',2,'LineWidth',2)

    theAlphaVal = 0.1;
    fvad = .8.*ones(size(theShade,1),1);
    fvad(length(pt_Tspan_smooth)) = theAlphaVal;
    fvad(end) = theAlphaVal;

    hold on
    patch(theShade(:,1),theShade(:,2),c1(1,:),'FaceAlpha',theAlphaVal,...
        'EdgeColor',c1(1,:),'EdgeAlpha','flat','LineStyle','-',...
        'FaceVertexAlphaData',fvad)
    %set(currplot,'Color',c1(1,:))
    hold off

    daTitle = sprintf('%3.0f nM',cellData.uconcs{expTime}(treatIter));
    title(daTitle)

    cAUC = sprintf('%3.0f',cellData.uconcs{expTime}(treatIter)*...
        cellData.ExposureTimes(expTime));
    dox_pss = sprintf('%3.1f\n',cellData.dox_pss{expTime}(treatIter));
    %text(LabelLocs_x(cntr), ax.YLim(2)*(LabelLocs_y(cntr)-.12), ['{\itAUC} = ' cAUC],'FontSize',fs)
    %text(LabelLocs_x(cntr), ax.YLim(2)*LabelLocs_y(cntr), ['{\itC_{B,max}} = ' dox_pss],'FontSize',fs)


    annLocs = allLocs(allLoc_loc(spi),:);
    an = annotation('textbox', annLocs,'String',...
        ['{\itC_{B,max}} = ' dox_pss '{\itAUC} = ' cAUC],...
        'FontSize',fs,'FitBoxToText','on');
    set(an,'FontName','Times')

    ax = gca;
    set(gca,'XLim',[-100 (cellData.Tspan{expTime}(end)+1*24)])
    set(gca,'YLim',[0 4e4])

    set(gca,'YTickLabel',sprintf('%2.1f\n',ax.YTick/1e4))
    set(gca,'FontSize',fs)
    text(-100, ax.YLim(2)*1.07, 'x10^4','FontSize',fs)
end
