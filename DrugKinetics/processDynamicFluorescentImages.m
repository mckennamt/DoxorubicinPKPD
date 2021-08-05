

function [outtie] = processDynamicFluorescentImages(Im_control,...
    Im_cell,cellMask)

%%read in 3 images
%Im_control;
%Im_cell;

% output corrected images


%Im_cell_current = double(imread(cell_flName))./maxInten;
            %Im_aif_current = double(imread(aif_flName))./maxInten;
            %
            %             validCount = 0;
            %             Im_cell = zeros(size(Im_cell_current));
            %             Im_aif = zeros(size(Im_aif_current));
            %
            %             % correct with histogram normalization
            %             for maf = allIm%allIm-1:allIm+1
            %
            %                 % aif and cell wells correction
            %                 flName = sprintf('Doxorubicin - n%06d.tif',maf);
            %                 %aif_flName = fullfile(controlFolder_AIF,flName);
            %                 aif_flName = fullfile(drugOnlyFolder,flName);
            %                 cell_flName = fullfile(cellFolder,flName);
            %
            %                 %Histogram correct consecutive images with drug
            %                 if (maf>=startDrug && maf<=endDrug) && ...
            %                         (allIm>=startDrug && allIm<=endDrug)
            %
            %                     validCount = validCount+1;
            %
            %                     Im_cell_maf = double(imread(cell_flName))./maxInten;
            %                     tpHist = hist(Im_cell_maf(:),davec);
            %                     Im_cell = Im_cell + histeq(Im_cell_current,tpHist);
            %
            %
            %                     Im_aif_maf = double(imread(aif_flName))./maxInten;
            %                     tpHist = hist(Im_aif_maf(:),davec);
            %                     Im_aif = Im_aif + histeq(Im_aif_current,tpHist);
            %
            %                 %Histogram correct consecutive images without drug
            %                 elseif (maf<startDrug || maf>endDrug) &&...
            %                         (allIm<startDrug || allIm>endDrug) &&...
            %                         (maf >=0)
            %
            %
            %                     validCount = validCount+1;
            %
            %                     Im_cell_maf = double(imread(cell_flName))./maxInten;
            %                     tpHist = hist(Im_cell_maf(:),davec);
            %                     Im_cell = Im_cell + histeq(Im_cell_current,tpHist);
            %
            %
            %                     Im_aif_maf = double(imread(aif_flName))./maxInten;
            %                     tpHist = hist(Im_aif_maf(:),davec);
            %                     Im_aif = Im_aif + histeq(Im_aif_current,tpHist);
            %
            %                 end
            %             end
            %
            %             Im_cell = Im_cell.*maxInten./validCount;
            %             Im_aif = Im_aif.*maxInten./validCount;
            %
            %              flNametl = sprintf('Transmitted Light - n%06d.tif',maf);
            %              cellTL_flName = fullfile(cellFolder,flNametl);
            %              Im_cellTL = double(imread(cellTL_flName));
            %
            %             % drug only wells
            %             flName = sprintf('Doxorubicin - n%06d.tif',allIm);
            %             drugOnly_flName = fullfile(drugOnlyFolder,flName);
            %             Im_control_current = double(imread(drugOnly_flName))./maxInten;
            %             Im_control = zeros(size(Im_control_current));
            %             validCount = 0;
            %             if allIm < startDrug
            %                 Im_control = Im_control_current.*maxInten;
            %             else
            %                 for maf = allIm%allIm-1:allIm+1
            %                     if (maf>=startDrug) && (allIm>=startDrug)
            %                         validCount = validCount + 1;
            %                         flName = sprintf('Doxorubicin - n%06d.tif',maf);
            %
            %                         Im_control_maf = double(imread(cell_flName))./maxInten;
            %                         tpHist = hist(Im_control_maf(:),davec);
            %                         Im_control = Im_control + histeq(Im_control_current,tpHist);
            %                     end
            %                 end
            %                 Im_control = Im_control.*maxInten./validCount;
            %             end