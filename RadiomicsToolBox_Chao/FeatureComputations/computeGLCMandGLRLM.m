function [GLCMnGLRLM] = computeTextureRadiomics(II, levels)
% levels, is replaced by Chao Huang according to 'prepareVolume.m' of 'mvallieres/radiomics'

 % m = floor(min(min(min(II))));
 % M = ceil(max(max(max(II))));
 % GLCM = computeGrayLevelCoocurrenceMatrix(II, m:M);
 GLCM = computeGrayLevelCoocurrenceMatrix(II, levels);
 GLCMRadiomics = computeGLCM(GLCM);
 % GLRLM = computeGrayLevelRunLengthMatrix(II, m:M);
 GLRLM = computeGrayLevelRunLengthMatrix(II, levels);
 GLRLMRadiomics = computeGLRLM(GLRLM);
 GLCMnGLRLM = catstruct(GLCMRadiomics,GLRLMRadiomics);
