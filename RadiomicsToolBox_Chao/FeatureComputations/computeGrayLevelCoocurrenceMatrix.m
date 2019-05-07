function glccm  = computeGrayLevelCoocurrenceMatrix(I, levels)

% this below part is actually uniformquantization. But I think it's reasonable to do both isotropic resampling and quantization before computing GLCM which is done in 'prepareVolume.m' of 'mvallieres/radiomics', in which some arguments could be customized. Thus, I prefer the later one for preprocess of the ROI.
% m = floor(min(min(min(I))));
% M = ceil(max(max(max(I))));    
% I = (I-m)/(M-m);
% I = I*10;
% I = round(I);

directions = [1 0 0; 1 1 0; -1 1 0; 0 1 0; 0 1 1; 0 -1 1; 0 0 1; 1 0 1; -1 0 1; -1 1 1; 1 -1 1; 1 1 1; -1 -1 1];
nDirections = 13;
nLevels = length(levels); % added by Chao: please find "levels" which equals [m:M], calculated in "computeTextureRadiomics.m".
glccm = zeros(nLevels, nLevels, nDirections);
nRow = size(I,1);
nCol = size(I,2);
nSli = size(I,3);
for k = 1:nSli
    for j = 1:nCol
        for i = 1:nRow
            currentPixelValue = I(i,j,k);
            rowIndex = find(levels == currentPixelValue);
            for l = 1:nDirections
                offset = [i+directions(l,1) j+directions(l,2) k+directions(l,3)];
                if (offset(1) < 1) || (offset(1) > nRow) || (offset(2) < 1) || (offset(2) > nCol) || (offset(3) < 1) || (offset(3) > nSli)
                    continue
                else
                    offsetPixelValue = I(offset(1), offset(2), offset(3));
                    colIndex = find(levels == offsetPixelValue);
                    glccm(rowIndex, colIndex, l) = glccm(rowIndex, colIndex, l)+1;
                end
            end
            
        end
    end
end