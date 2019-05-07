function [Radiomics] = computeShapeSize(ROIonly,pixelW,sliceS)
% The Volume and Surface area, etc were computed based on algorithms defined by Hugo Aerts 2014 Nature Communication Paper ('Decoding...'). Some codes were adapted from 'computeImageRadiomics.m' of 'chihyanghsu0805/Radiomics' and 'mvallieres/radiomics' folder.
% number of shape and size features (three-dimensional features): 9. updated Dec16, 2017
% N = 10

[x,y,z] = find3d(ROIonly); % returns the i,j,k indices of non-zero entries. adapted from
% 'computeBoundingBox.m' of <https://github.com/mvallieres/radiomics/>
pixelList = [x',y',z'];

% compute volume
V = sum(~isnan(ROIonly(:)))*pixelW*pixelW*sliceS;
maskBox = double(~isnan(ROIonly)); % get 1 where ROIonly ~= NaN, 0 where ROIonly = NaN.

% compute surface area, refer: https://stackoverflow.com/questions/25747918/how-to-calculate-the-area-by-triangulating-a-3d-object-in-matlab
FV = isosurface(maskBox, 0.5); % get an isosurface where values=0.5
vertices = FV.vertices; % shape: num_vertices,3. 3 means x,y,z.
faces = FV.faces; % shape: num_faces, 3. 3 means A,B,C vertices indices (range: 1 to num_vertices) of a triangle.
scale = pixelW*pixelW*sliceS;
a = vertices(faces(:, 2), :) - vertices(faces(:, 1), :); % an edge vector of triangle
b = vertices(faces(:, 3), :) - vertices(faces(:, 1), :); % an edge vector of triangle
cro = cross(a, b, 2);
scaled_cro = cro*scale*scale; % cro was computed assuming the pixel resolution is 1*1*1;
A = 1/2*sum(sqrt(sum(scaled_cro.^2, 2)));

% compute max diameter
surf_pixels_cor = [vertices(:,1) * pixelW, vertices(:,2) * pixelW, vertices(:,3) * sliceS]; % iscolumn(vertices(:,1))==1

 % compute maximum diameter based on Hugo Aerts Nat Commu paper
step = 30000;% 50000;
Num = ceil(length(surf_pixels_cor)/step);
if Num ==1
    Radiomics.shapeSize_Maximum3DDiameter = max(pdist(surf_pixels_cor, 'euclidean'));
else
    % to handle with too many points in which case PC might be out of memory 
    mD = 0;
    N_prod = carte_prod(1:Num, 1:Num);
    parfor k=1:length(N_prod)
        num = N_prod(k,:);
        i = num(1); j = num(2);
        start1 = (i-1)*step+1;
        end1 = min(i*step, length(surf_pixels_cor));
        start2 = (j-1)*step+1;
        end2 = min(j*step, length(surf_pixels_cor));
        pd = pdist2(surf_pixels_cor(start1:end1,:), surf_pixels_cor(start2:end2, :));
        mD = max(mD, max(pd(:)));
    end
    Radiomics.shapeSize_Maximum3DDiameter = mD;
end


% Algorithms same to Aerts' paper
Radiomics.shapeSize_Compactness1 = V/sqrt(pi)/A.^(2/3);
Radiomics.shapeSize_Compactness2 = 36*pi*V.^2/A.^3;
R = (V*3/4/pi).^(1/3);
Radiomics.shapeSize_SphericalDisproportion = A/4/pi/R.^2;
Radiomics.shapeSize_Sphericity = pi.^(1/3)*(6*V).^(2/3)/A;
Radiomics.shapeSize_Volume = V; 
Radiomics.shapeSize_SurfaceArea = A;
Radiomics.shapeSize_SurfaceVolumeRatio = A/V;

% Eccentricity and Solidity are from <https://github.com/mvallieres/radiomics/>
if length(size(ROIonly)) ==2
Radiomics.shapeSize_Eccentricity = 0;
Radiomics.shapeSize_Solidity=1;
else
Radiomics.shapeSize_Eccentricity = getEccentricity(ROIonly,pixelW,sliceS); 
Radiomics.shapeSize_Solidity = getSolidity(ROIonly,pixelW,sliceS);
end

