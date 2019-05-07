function [wavelet] = computeWavelet(volume,mask,scanType,pixelW,sliceS,R,scale,textType,quantAlgo,Ng) 
if strcmp(quantAlgo,'Lloyd')
    quantization = @(x,y) lloydQuantization(x,y);
elseif strcmp(quantAlgo,'Equal')
    quantization = @(x,y) equalQuantization(x,y);
elseif strcmp(quantAlgo,'Uniform')
    quantization = @(x,y) uniformQuantization(x,y);
else
    error('Error with quantization algorithm input. Must either be ''Equal'' or ''Lloyd'' or ''Uniform''')
end

%% non-texture features
R = 1; % states no ROI fill and no waveletfilt as in prepareVolume().
[I,levels,ROIbox,maskBox] = prepareVolume(volume,mask,scanType,pixelW,sliceS,R,scale,'No_preprocessing');

if sum(isnan(I(:)))
    I = fillBox(I); % Necessary in cases we have a ROI box containing NaN's.
end
    
waveDec = wavedec3(I, 1, 'coif1');

% the order is stated in help-wavedec3.
% waveDec.dec{1}; %% LLL
% waveDec.dec{2}; %% HLL
% waveDec.dec{3}; %% LHL
% waveDec.dec{4}; %% HHL
% waveDec.dec{5}; %% LLH
% waveDec.dec{6}; %% HLH
% waveDec.dec{7}; %% LHH
% waveDec.dec{8}; %% HHH
waveletNames = {'LLL','HLL','LHL','HHL','LLH','HLH','LHH','HHH'};

firstOrderWavelet = struct;

for i = 1:8
    waveletName = char(waveletNames(i));
    II = waveDec.dec{i};
    % firstOrder/intensity-based/global features.
    firstOrder = computeFirstOrder(II);
    firstOrder = renameStructFields(firstOrder, [waveletName,'_'], 'wavelet');
    firstOrderWavelet = catstruct(firstOrderWavelet, firstOrder);
end

%% texture features
R = 1; % states no ROI fill and no waveletfilt as in prepareVolume().
[I,levels,ROIbox,maskBox] = prepareVolume(volume,mask,scanType,pixelW,sliceS,R,scale,'No_preprocessing');

if sum(isnan(I(:)))
    I = fillBox(I); % Necessary in cases we have a ROI box containing NaN's.
end
    
waveDec = wavedec3(I, 1, 'coif1');

% the order is stated in help-wavedec3.
% waveDec.dec{1}; %% LLL
% waveDec.dec{2}; %% HLL
% waveDec.dec{3}; %% LHL
% waveDec.dec{4}; %% HHL
% waveDec.dec{5}; %% LLH
% waveDec.dec{6}; %% HLH
% waveDec.dec{7}; %% LHH
% waveDec.dec{8}; %% HHH
waveletNames = {'LLL','HLL','LHL','HHL','LLH','HLH','LHH','HHH'};

GLCMnGLRLMWavelet = struct;

for i = 1:8
    waveletName = char(waveletNames(i));
    II = waveDec.dec{i};
    
    %%%%
    ROIbox = II;

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
            ROIbox = imresize3D(ROIbox,[],[round(double(size(ROIbox,1))*a),round(double(size(ROIbox,2))*b),round(double(size(ROIbox,3))*c)],'cubic','fill');
        end
    elseif numel(size(ROIbox))==2
        if a + b ~= 2 % If false, no resampling is needed
            ROIbox = imresize(ROIbox,[round(double(size(ROIbox,1))*a),round(double(size(ROIbox,2))*b)],'cubic','Antialiasing',true);
        end
    end
    
    ROIonly = ROIbox;
    [ROIonly,levels] = quantization(ROIonly,Ng);
    
    %%%%
    II_matrix = ROIonly;
    
    GLCM = computeGrayLevelCoocurrenceMatrix(II_matrix, levels);
    GLCMRadiomics = computeGLCM(GLCM);
    
    GLRLM = computeGrayLevelRunLengthMatrix(II_matrix, levels);
    GLRLMRadiomics = computeGLRLM(GLRLM);
    
    GLCMnGLRLM = catstruct(GLCMRadiomics,GLRLMRadiomics);
    GLCMnGLRLM = renameStructFields(GLCMnGLRLM, [waveletName,'_'], 'wavelet');
    GLCMnGLRLMWavelet = catstruct(GLCMnGLRLMWavelet, GLCMnGLRLM);
end

%% collect all features
    wavelet = catstruct(firstOrderWavelet, GLCMnGLRLMWavelet);
end
