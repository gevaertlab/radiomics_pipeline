function [CoLlAGe_theta, CoLlAGe_fai] = computeCoLlAGe(volume,mask,scanType,pixelW,sliceS,scale)
% by Chao. Updated: Oct27, 2018.
% Compute CoLlAGe features, based on:
% Paper:
%     Prasanna et al. 2016, Co-occurrence of Local Anisotropic
%     Gradient Orientations (CoLlAGe): A
%     new radiomics descriptor, Nature Scientific Reports.
% & some code refer to: 
%     https://github.com/Svdvoort/PREDICTFastr/

%Args:
%    volume: not preprocessed volume. Will be properly preprocessed in this
%    function.

% the following are hyperparameters which need to be optimized. Here, for nbins and window_size, we
% list more than potential values listed in the paper.

nbins = [10 20 30]; % symbol: v, the entropy histogram is divided into bin size v. 
window_size = [5 7 9]; % symbol: N, a local neighborhood of N voxels defined in a window
omega = 64; % symbol: w, discretization factor for entropy map
% wsize = 1; dist = 2; % To get GLCM matrix,instead of skimage.feature.greycomatrix in python, we used the
% computeGrayLevelCoocurrenceMatrix() which applies to every pair of levels(i) and levels(j), and
% totally 13 directions on x,y,z. Details go to computeGrayLevelCoocurrenceMatrix()

% COMPUTATION OF THE SMALLEST BOX CONTAINING THE ROI
[boxBound] = computeBoundingBox(mask);
maskBox = mask(boxBound(1,1):boxBound(1,2),boxBound(2,1):boxBound(2,2),boxBound(3,1):boxBound(3,2));
ROIbox = volume(boxBound(1,1):boxBound(1,2),boxBound(2,1):boxBound(2,2),boxBound(3,1):boxBound(3,2));

% PRE-PROCESSING OF ROI BOX
ROIbox = double(ROIbox);
if strcmp(scanType,'PETscan') || strcmp(scanType,'PTscan')
    ROIbox = sqrt(ROIbox);
elseif strcmp(scanType,'MRscan')
    ROIonly = ROIbox;
    ROIonly(~maskBox) = NaN;
    temp = CollewetNorm(ROIonly);
    maskBox(isnan(temp)) = 0;
end

% fill NaNs
if sum(isnan(ROIbox(:)))
    ROIbox = fillBox(ROIbox); % Necessary in cases we have a ROI box containing NaN's.
end

% ISOTROPIC RESAMPLING
flagPW = 0;
if strcmp(scale,'pixelW')
    flagPW = 1;
end
if flagPW
    a = 1;
    b = 1;
    c = sliceS/pixelW;
else
    a = pixelW/scale;
    b = pixelW/scale;
    c = sliceS/scale;
end
if numel(size(ROIbox))==3
    if a + b + c ~= 3 % If false, no resampling is needed
        maskBox = imresize3D(maskBox,[],[round(double(size(maskBox,1))*a),round(double(size(maskBox,2))*b),round(double(size(maskBox,3))*c)],'nearest','fill');
        ROIbox = imresize3D(ROIbox,[],[round(double(size(ROIbox,1))*a),round(double(size(ROIbox,2))*b),round(double(size(ROIbox,3))*c)],'cubic','fill');
    end
elseif numel(size(ROIbox))==2
    if a + b ~= 2 % If false, no resampling is needed
        maskBox = imresize(maskBox,[round(double(size(maskBox,1))*a),round(double(size(maskBox,2))*b)],'nearest');
        ROIbox = imresize(ROIbox,[round(double(size(ROIbox,1))*a),round(double(size(ROIbox,2))*b)],'cubic','Antialiasing',true);
    end
end

% compute gradients
[gradient_x, gradient_y, gradient_z] = imgradientxyz(ROIbox); % ROIbox will not have NaN, and should be able to compute gradients for every pixel of interest

% determine ROIonly
ROIonly = ROIbox;
ROIonly(~maskBox) = NaN;
ROIonly(maskBox<0) = NaN; 
[x_ids, y_ids, z_ids] = find3d(ROIonly); % returns the i,j,k indices of non-zero entries.
xmin = min(x_ids(:)); xmax = max(x_ids(:));
ymin = min(y_ids(:)); ymax = max(y_ids(:));
zmin = min(z_ids(:)); zmax = max(z_ids(:));
pixel_num = (xmax-xmin+1)*(ymax-ymin+1)*(zmax-zmin+1);


% initiates
E = struct; % E, e.g. entropy
for v = nbins
    E.(strcat('v',num2str(v))) = struct;
end


for v = nbins
    for N = window_size
        theta_E.(strcat('v',num2str(v))).(strcat('N',num2str(N))) = zeros((xmax-xmin+1)*(ymax-ymin+1)*(zmax-zmin+1),1);
        fai_E.(strcat('v',num2str(v))).(strcat('N',num2str(N))) = zeros((xmax-xmin+1)*(ymax-ymin+1)*(zmax-zmin+1),1);
        
        theta = zeros(xmax-xmin+1, ymax-ymin+1, zmax-zmin+1); % init dominant orientations across(X,Y).
        fai = zeros(xmax-xmin+1, ymax-ymin+1, zmax-zmin+1); % init dominant orientations across(X,Y,Z).
        % vectorize to apply parfo-loop
        theta_tmp = zeros(pixel_num,1);
        fai_tmp = zeros(pixel_num,1);
        
        c_prod = carte_prod(xmin:xmax,ymin:ymax,zmin:zmax);
        fprintf('CoLlAGe: computing theta and fai...v:%d, N:%d -----', v, N)
        tic % start to cal time
        parfor i = 1:length(c_prod)
            c = c_prod(i,:,:);
            x = c(1); y = c(2); z = c(3);
            
            % Construct gradient vector for every pixel per cell
            coord_prod = carte_prod((x-floor(N/2)):(x+floor(N/2)), (y-floor(N/2)):(y+floor(N/2)), (z-floor(N/2)):(z+floor(N/2)));
            Fx = zeros(N^3,1); Fy = zeros(N^3,1); Fz = zeros(N^3,1);
            for j = 1:length(coord_prod)
                coord_N = coord_prod(j,:,:);
                % NOTE: Not covered in paper how to treat boundary pixels!
                % If out of bounds, use boundary pixels
                xc = min(max(coord_N(1), xmin),xmax);
                yc = min(max(coord_N(2), ymin),ymax);
                zc = min(max(coord_N(3), zmin),zmax);
%                 fprintf('xc %f, yc %f, zc %f \n', [xc yc zc]);
                % Append gradients
                Fx(j,1) = gradient_x(xc, yc, zc);
                Fy(j,1) = gradient_y(xc, yc, zc);
                Fz(j,1) = gradient_z(xc, yc, zc);
            end
            % Store gradients in single matrix
            F = zeros(N^3, 3);
            F(:,1) = Fx; F(:,2) = Fy; F(:,3) = Fz;
            
            % Get dominant gradient through SVD
            [U, s, V] = svd(F); % V shape: (3,3).?(ck11), ?(ck12), and ?(ck13) of svd will be extracted later
%             fprintf('x-xmin %f, y-ymin %f, z-zmin %f \n', [x-xmin, y-ymin, z-zmin]);
%             theta(x-xmin, y-ymin, z-zmin) = atan(V(2,1)/V(1,1));
%             fai(x-xmin, y-ymin, z-zmin) = atan(V(3,1)/sqrt((V(2,1)^2+V(1,1)^2)));
%             theta(x,y,z) = atan(V(2,1)/V(1,1));
%             fai(x,y,z) = atan(V(3,1)/sqrt((V(2,1)^2+V(1,1)^2)));
            theta_tmp(i) = atan(V(2,1)/V(1,1));
            fai_tmp(i) = atan(V(3,1)/sqrt((V(2,1)^2+V(1,1)^2)));
        end
        toc % end to cal time
        theta = reshape(theta_tmp, [xmax-xmin+1, ymax-ymin+1, zmax-zmin+1]);
        fai = reshape(fai_tmp, [xmax-xmin+1, ymax-ymin+1, zmax-zmin+1]);
        % Discretize theta and fai analogue to parsanna et al.
        theta = omega*(theta-min(theta(:)))/(max(theta(:))-min(theta(:)));
        theta = uint8(theta);
        fai = omega*(fai-min(fai(:)))/(max(fai(:))-min(fai(:)));
        fai = uint8(fai);
        
        % entropy. Again, computations for each pixel
        theta_E_tmp = zeros(xmax-xmin+1, ymax-ymin+1, zmax-zmin+1);
        fai_E_tmp = zeros(xmax-xmin+1, ymax-ymin+1, zmax-zmin+1);
        % vectorize to apply parfo-loop
        theta_E_tmp_tmp = zeros(pixel_num,1);
        fai_E_tmp_tmp = zeros(pixel_num,1);
        
        c_prod = carte_prod(1:size(theta,1), 1:size(theta,2), 1:size(theta,3));
        fprintf('CoLlAGe: computing entropy...v:%d, N:%d -----', v, N)
        tic % start to cal time
        parfor i=1:length(c_prod)
            c = c_prod(i,:,:);
            x = c(1); y = c(2); z = c(3);
            
            % If out of bounds of mask, use boundary pixels
            xstart = max(x-floor(N/2), 1); xend = min(x+floor(N/2), size(theta,1));
            ystart = max(y-floor(N/2), 1); yend = min(y+floor(N/2), size(theta,2));
            zstart = max(z-floor(N/2), 1); zend = min(z+floor(N/2), size(theta,3));
            
            % Compute gray level co-occurence matrix
            theta_window = theta(xstart:xend, ystart:yend, zstart:zend);
            fai_window = fai(xstart:xend, ystart:yend, zstart:zend);
            theta_GLCM_matrix = computeGrayLevelCoocurrenceMatrix(theta_window, omega); % GLCM shape: zeros(nLevels, nLevels, nDirections);
            fai_GLCM_matrix = computeGrayLevelCoocurrenceMatrix(fai_window, omega);
            
            % Compute entropy
            theta_e = entropy(theta_GLCM_matrix);
            fai_e = entropy(fai_GLCM_matrix);
            if theta_e<0 && isinf(theta_e)
                theta_e = 1e-12;
            end
            if fai_e<0 && isinf(fai_e)
                fai_e = 1e-12;
            end
            theta_E_tmp_tmp(i) = theta_e;
            fai_E_tmp_tmp(i) = fai_e;
%             theta_E_tmp(x,y,z) = theta_e;
%             fai_E_tmp(x,y,z) = fai_e;
        end
        toc % end to cal time
        theta_E_tmp = reshape(theta_E_tmp_tmp, [xmax-xmin+1, ymax-ymin+1, zmax-zmin+1]);
        fai_E_tmp = reshape(fai_E_tmp_tmp, [xmax-xmin+1, ymax-ymin+1, zmax-zmin+1]);
        theta_E.(strcat('v',num2str(v))).(strcat('N',num2str(N)))(:,1) = theta_E_tmp(:);
        fai_E.(strcat('v',num2str(v))).(strcat('N',num2str(N)))(:,1) = fai_E_tmp(:); 
    end
end

% compute histogram
fprintf('CoLlAGe: histgraming entropy...')
tic
for v = nbins
    theta_E_nbins = theta_E.(strcat('v',num2str(v)));
    fai_E_nbins = fai_E.(strcat('v',num2str(v)));
    for N = window_size
        
        theta_E_ws = theta_E_nbins.(strcat('N',num2str(N))); % E_ws for every window size
        fai_E_ws = fai_E_nbins.(strcat('N',num2str(N)));
        theta_hist = histogram(theta_E_ws, uint8(v));
        theta_hist_values = theta_hist.Values;
        fai_hist = histogram(fai_E_ws, uint8(v));
        fai_hist_values = fai_hist.Values;
        
        for bin=1:v
            CoLlAGe_theta.(strcat('CoLlAGe_thetav',num2str(v),'N',num2str(N), 'bin',num2str(bin))) = theta_hist_values(bin);
            CoLlAGe_fai.(strcat('CoLlAGe_faiv',num2str(v),'N',num2str(N), 'bin',num2str(bin))) = fai_hist_values(bin);
        end
    end
end
toc


end


