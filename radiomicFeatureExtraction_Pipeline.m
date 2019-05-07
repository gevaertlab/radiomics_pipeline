function [radioFeatures] = radiomicFeatureExtraction_Pipeline(DataVolume,MaskVolume,Para)
%% UNTITLED Summary of this function goes here
% xxxx features in total.
% updated: Apr 4th, 2019.

%% set up
% paras to pass to perpareVolume()
volume = DataVolume;
mask = MaskVolume;
scanType = Para.type;
pixelW = Para.pixelW;
sliceS = Para.sliceS;
R = Para.R;
scale = Para.scale;
quantAlgo = Para.quantAlgo;
Ng = Para.Ng;

%% Prepare volume to compute intensity-based features.
ROI_only_no_preproc = prepareVolume(volume,mask,scanType,pixelW,sliceS,R,scale,'No_preprocessing');

%% Group#1: first order (N=21)
fprintf('Extracting first order... \n');
[firstOrder] = computeFirstOrder(ROI_only_no_preproc);

%% Group#2: shape&size (N=10)
fprintf('Extracting shape&size features... \n');
[shapeFeatures] = computeShapeSize(ROI_only_no_preproc, pixelW,sliceS);

%% Group#3: Global histogram features (N = 17)
fprintf('Extracting Global histogram features... \n');
Nbins=50;
GlobalHistogramFeatures = getGlobalHistogram(ROI_only_no_preproc,Nbins);


%% Group#4: texture
% fprintf('\nExtracting texture features... \n');

%% subgroup1: matrix based (total N=52). GLCM(N=23), GLRLM(N=11), GLSZM(N=13), NGTDM(N=5).
fprintf('Extracting texture: matrix-based features... \n');
% To prepare volume for matrix-based texture features extraction.
% m stands for 'Matrix'
[ROIonly_m,levels_m,ROIbox_m,maskBox_m] = prepareVolume(volume,mask,scanType,pixelW,sliceS,R,scale,'Matrix',quantAlgo,Ng); % from 'mvallieres/radiomics' folder.

% subgroup3.1: matrix-based texture features % n approximate workaround for dealing with
% non-rectangular ROI's would be setting the pixels outside the ROI to 0 and using the funtion haralick from mahotas library.
% from "chihyanghsu0805/Radiomics", including quantization before computing GLCM
% matrix and GLRLM matrix. Not using methods from 'mvallieres/radiomics' to compute GLCM and GLRLM because they are different from Hugo Aerts' Nat Com description.
GLCM = computeGrayLevelCoocurrenceMatrix(ROIonly_m, levels_m);
GLCMRadiomics = computeGLCM(GLCM);

GLRLM = computeGrayLevelRunLengthMatrix(ROIonly_m, levels_m);
GLRLMRadiomics = computeGLRLM(GLRLM);

[GLSZM_matrix] = getGLSZM(ROIonly_m,levels_m); % 'mvallieres/radiomics'
[GLSZM_tmp] = getGLSZMtextures(GLSZM_matrix);
[GLSZM] = renameStructFields(GLSZM_tmp,'_','GLSZM');

[NGTDM_matrix,countValid] = getNGTDM(ROIonly_m,levels_m); % 'mvallieres/radiomics'
[NGTDM_tmp] = getNGTDMtextures(NGTDM_matrix,countValid);
[NGTDM] = renameStructFields(NGTDM_tmp,'_','NGTDM');

matrixTextures = catstruct(GLCMRadiomics,GLRLMRadiomics, GLSZM, NGTDM);

%% subgroup2: CoLlAGe (N=540)
% fprintf('Extracting texture: CoLlAGe features... \n');
% [CoLlAGe_theta, CoLlAGe_fai] = computeCoLlAGe(volume,mask,scanType,pixelW,sliceS,scale);

%% Group#5:  filter based
%% subgroup1: wavelet filter transformed features (N=464)
fprintf('Extracting wavelet features... \n');
% Only wavelet features for first order, GLCM and GLRLM now. If necessary,
% add codes for GLSZM, NGTDM
[wavelet] = computeWavelet(volume,mask,scanType,pixelW,sliceS,R,scale,'',quantAlgo,Ng);

% return all the features
radioFeatures = catstruct(firstOrder, shapeFeatures, GlobalHistogramFeatures, matrixTextures, wavelet);
% radioFeatures = catstruct(firstOrder, shapeFeatures, GlobalHistogramFeatures, matrixTextures, HOGfeaturesFinal, LBPfeaturesFinal, wavelet, GaborFeatures);

end % end of function definition

