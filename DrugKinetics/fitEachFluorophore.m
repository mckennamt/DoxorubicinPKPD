function unmixed = fitEachFluorophore(measured, T)

% unmixed: image size x number fluorophores
% measured: image size x number channels
% T: matrix relating true values to measured values (number fluorophores x
% number image channels)

% estimate true image for each fluorophore
% most likely from channel with highest signal
estimatedImages = zeros(size(measured,1), size(measured,2), size(T,1));
for fluorIter = 1:size(T,1)
    estimatedImages(:,:,fluorIter) = measured(:,:,T(fluorIter,:) == 1);
end

unmixed = zeros(size(estimatedImages));

% each image too big to simultaneously fit
% need to separate into separate images (10x10)
% need a for loop to do this
kernalSize = 16; %square
halfKernal = kernalSize/2;

kernalCenterRow = kernalSize/2:(kernalSize):imSize(1);
kernalCenterCol = kernalSize/2:(kernalSize):imSize(2);

for rowIter = 1:length(kernalCenterRow)
    for colIter = 1:length(kernalCenterCol)
        
        cr = kernalCenterRow(rowIter);
        cc = kernalCenterCol(colIter);
        
        % get image sub-sections
        estimatedImages_ss = estimatedImages(cr-(halfKernal-1):cr+halfKernal,...
            cc-(halfKernal-1):cc+halfKernal,:);
        
        measured_ss = measured(cr-(halfKernal-1):cr+halfKernal,...
            cc-(halfKernal-1):cc+halfKernal,:);
        
        % reshape images into arrays
        % num pixels x num fluorophores
        c_estimate = reshape(estimatedImages_ss,kernalSize^2,size(T,1));
        
        % num pixels x num channels
        c_measured = reshape(measured_ss,kernalSize^2,size(T,2));
        
        Options=optimset('Display','none','TolFun',1e-15,'TolX',1e-5);
        c_optimized = lsqnonlin(@(x) ...
            unmixImages(x,c_measured,T), c_estimate,...
            [],[],Options);
        
        unmixed(cr-(halfKernal-1):cr+halfKernal,...
            cc-(halfKernal-1):cc+halfKernal,:) = ...
            reshape(c_optimized,kernalSize,kernalSize,size(T,1));
        
    end
end

end


function Y = unmixImages(estimatedImages, measured, T)

est_measured = estimatedImages*T;

% ||Ymeas - Yest||2 + reg*||smoothness||2
Y = measured(:) - est_measured(:);

end