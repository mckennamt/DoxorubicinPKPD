
% Load bootstrap data
allCellLine = {'SUM149','MDAMB468', 'MDAMB231'};
CellLine = allCellLine{3};
load(fullfile(pwd, ['BootstrapFits_' CellLine '.mat']))
%%
c1 = distinguishable_colors(1);

% iterate through each exposure time for each cell line
for expTimes = 1:size(bestFitModels,3)
    
    % iter through each concentration
    for treatIter = 1:size(bestFitModels,1)
        
        cFits = cat(2, bestFitModels_smooth{treatIter, :, expTimes});
        
        allTspan = cellData.Tspan{expTimes};
        postTreat_TP = cellData.DrugAdded_Tp + 1;
        pt_Tspan = allTspan(postTreat_TP:end);
        pt_Tspan_smooth = linspace(pt_Tspan(1), pt_Tspan(end), size(cFits,1));
        
        inds_95CI = ceil(.025*numBootstrap);
        theShade = zeros(length(pt_Tspan_smooth)*2,2);
        theShade(:,1) = [pt_Tspan_smooth'; pt_Tspan_smooth(end:-1:1)'];
        for makeShade = 1:length(pt_Tspan_smooth)
            srted = sort(cFits(makeShade,:));
            theShade(makeShade,2) = srted(inds_95CI);
        end
        
        inds_95CI = floor(.975*numBootstrap);
        shadeIndex = length(pt_Tspan_smooth);
        for makeShade = length(pt_Tspan_smooth):-1:1
            shadeIndex = shadeIndex+1;
            srted = sort(cFits(makeShade,:));
            theShade(shadeIndex,2) = srted(inds_95CI);
        end
        
        
        cellCounts = cellData.allYall{expTimes}{treatIter+1};
        daStd = std(cellCounts');
        daMean = mean(cellCounts');
        nSamp = size(cellCounts,2);
        
        bnds = tinv([0.025, 0.975],nSamp-1)' * daStd/sqrt(nSamp);
        ub = bnds(2,:);
        lb = -bnds(1,:);
        
        figure(expTimes);clf;
        errorbar(allTspan,daMean,...
            lb,ub,...
            'ks','MarkerSize',2,'LineWidth',2)
        
        hold on
        theAlphaVal = .1;
        patch(theShade(:,1),theShade(:,2),c1(1,:),'FaceAlpha',theAlphaVal,...
            'EdgeColor',c1(1,:),'EdgeAlpha',.2,'LineStyle','-')
        %set(currplot,'Color',c1(1,:))
        hold off
        
        
        pause
        
    end
    
end