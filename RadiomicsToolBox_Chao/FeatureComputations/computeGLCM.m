function Radiomics = computeGLCM(glcm)
% glcm: glcm matrix
nMatrices = size(glcm,3);
nLevels = size(glcm,2);

Radiomics.GLCM_Autocorrelation = zeros(nMatrices,1);
Radiomics.GLCM_ClusterProminence = zeros(nMatrices,1);
Radiomics.GLCM_ClusterShade = zeros(nMatrices,1);
Radiomics.GLCM_ClusterTendency = zeros(nMatrices,1);
Radiomics.GLCM_Contrast = zeros(nMatrices,1);
Radiomics.GLCM_Correlation = zeros(nMatrices,1);
Radiomics.GLCM_DifferenceEntropy = zeros(nMatrices,1);
Radiomics.GLCM_DifferenceVariance = zeros(nMatrices,1); % not in Aerts 2014, but in IBSI by Zwanenburg et al.
Radiomics.GLCM_Dissimilarity = zeros(nMatrices,1);
Radiomics.GLCM_Energy = zeros(nMatrices,1);
Radiomics.GLCM_Entropy = zeros(nMatrices,1);
Radiomics.GLCM_Homogeneity1 = zeros(nMatrices,1);
Radiomics.GLCM_Homogeneity2 = zeros(nMatrices,1);
Radiomics.GLCM_IMC1 = zeros(nMatrices,1);
Radiomics.GLCM_IMC2 = zeros(nMatrices,1);
Radiomics.GLCM_IDN = zeros(nMatrices,1);
Radiomics.GLCM_IDMN = zeros(nMatrices,1);
Radiomics.GLCM_InverseVariance = zeros(nMatrices,1);
Radiomics.GLCM_MaximumProbability = zeros(nMatrices,1);
Radiomics.GLCM_SumAverage = zeros(nMatrices,1);
Radiomics.GLCM_SumEntropy = zeros(nMatrices,1);
Radiomics.GLCM_SumVariance = zeros(nMatrices,1);
Radiomics.GLCM_Variance = zeros(nMatrices,1);

for k = 1:nMatrices
    P = glcm(:,:,k)/sum(sum(glcm(:,:,k)));
    mu = mean2(P);
    px = sum(P,1);
    py = sum(P,2);
    mux = sum((1:nLevels)*P);
    muy = sum((1:nLevels)*P');
    stdx = sum(((1:nLevels)-mux).^2*P);
    stdy = sum(((1:nLevels)-muy).^2*P');
    HX = -(px+eps)*log(px+eps)'; %HX = entropy(px);
    HY = -(py+eps)'*log(py+eps); %HX = entropy(px);
    temp = reshape(P,nLevels*nLevels,1);
    H = -(temp+eps)'*log(temp+eps);    
    HXY1 = 0;
    HXY2 = 0;
    pxminusy = zeros(nLevels,1);
    pxplusy = zeros(2*nLevels-1,1);
    for i = 1:nLevels
        for j = 1:nLevels
            Radiomics.GLCM_Autocorrelation(k,1) = i*j*P(i,j)+Radiomics.GLCM_Autocorrelation(k,1);
            Radiomics.GLCM_ClusterProminence(k,1) = (i+j-mux-muy)^4*P(i,j)+Radiomics.GLCM_ClusterProminence(k,1);
            Radiomics.GLCM_ClusterShade(k,1) = (i+j-mux-muy)^3*P(i,j)+Radiomics.GLCM_ClusterShade(k,1);
            Radiomics.GLCM_ClusterTendency(k,1) = (i+j-mux-muy)^2*P(i,j)+Radiomics.GLCM_ClusterTendency(k,1);
            Radiomics.GLCM_Contrast(k,1) = (i-j)*(i-j)*P(i,j)+Radiomics.GLCM_Contrast(k,1);
            if stdx*stdy == 0
                Radiomics.GLCM_Correlation(k,1) = 0+Radiomics.GLCM_Correlation(k,1);
            else
                Radiomics.GLCM_Correlation(k,1) = ((i-mux)*(j-muy)*P(i,j))/sqrt(stdx*stdy)+Radiomics.GLCM_Correlation(k,1); % Different Definition
            end
            Radiomics.GLCM_Dissimilarity(k,1) = abs(i-j)*P(i,j)+Radiomics.GLCM_Dissimilarity(k,1);
            Radiomics.GLCM_Homogeneity1(k,1) = P(i,j)/(1+abs(i-j))+Radiomics.GLCM_Homogeneity1(k,1);
            Radiomics.GLCM_Homogeneity2(k,1) = P(i,j)/(1+abs(i-j)^2)+Radiomics.GLCM_Homogeneity2(k,1);
            Radiomics.GLCM_IDMN(k,1) = P(i,j)/(1+(abs(i-j)^2/nLevels^2))+Radiomics.GLCM_IDMN(k,1);
            Radiomics.GLCM_IDN(k,1) = P(i,j)/(1+(abs(i-j)/nLevels))+Radiomics.GLCM_IDN(k,1);
            Radiomics.GLCM_Variance(k,1) = (i-mu)^2*P(i,j)+Radiomics.GLCM_Variance(k,1);
            
            if i == j
            else
                Radiomics.GLCM_InverseVariance(k,1) = P(i,j)/abs(i-j)^2+Radiomics.GLCM_InverseVariance(k,1);
            end
            HXY1 = HXY1+P(i,j)*log(px(i)*py(j)+eps);
            HXY2 = HXY2+px(i)*py(j)*log(px(i)*py(j)+eps);
            
            kk = abs(i-j)+1;
            ll = i+j-1;
            pxminusy(kk,1) = pxminusy(kk,1)+P(i,j);
            pxplusy(ll,1) = pxplusy(ll,1)+P(i,j);
        end
    end
    Radiomics.GLCM_Energy(k,1) = sum(sum(P.^2));
    Radiomics.GLCM_Entropy(k,1) = H;
    Radiomics.GLCM_DifferenceEntropy(k,1) = -(pxminusy+eps)'*log(pxminusy+eps);
    Radiomics.GLCM_DifferenceVariance(k,1) = (0:nLevels-1).^2*pxminusy;
    HXY1 = -HXY1;
    HXY2 = -HXY2;
    Radiomics.GLCM_IMC1(k,1) = (H-HXY1)/max(HX,HY);
    Radiomics.GLCM_IMC2(k,1) = sqrt(1-exp(-2*(HXY2-H)));
    Radiomics.GLCM_MaximumProbability(k,1) = max(max(P));
    Radiomics.GLCM_SumAverage(k,1) = (2:2*nLevels)*pxplusy;
    Radiomics.GLCM_SumEntropy(k,1) = -(pxplusy'+eps)*log(pxplusy+eps);
    Radiomics.GLCM_SumVariance(k,1) = ((2:2*nLevels)-Radiomics.GLCM_SumEntropy(k,1)).^2*pxplusy;
end

Radiomics.GLCM_Autocorrelation = mean(Radiomics.GLCM_Autocorrelation);
Radiomics.GLCM_ClusterProminence = mean(Radiomics.GLCM_ClusterProminence);
Radiomics.GLCM_ClusterShade = mean(Radiomics.GLCM_ClusterShade);
Radiomics.GLCM_ClusterTendency = mean(Radiomics.GLCM_ClusterTendency);
Radiomics.GLCM_Contrast = mean(Radiomics.GLCM_Contrast);
Radiomics.GLCM_Correlation = mean(Radiomics.GLCM_Correlation);
Radiomics.GLCM_DifferenceEntropy = mean(Radiomics.GLCM_DifferenceEntropy);
Radiomics.GLCM_DifferenceVariance = mean(Radiomics.GLCM_DifferenceVariance);
Radiomics.GLCM_Dissimilarity = mean(Radiomics.GLCM_Dissimilarity);
Radiomics.GLCM_Energy = mean(Radiomics.GLCM_Energy);
Radiomics.GLCM_Entropy = mean(Radiomics.GLCM_Entropy);
Radiomics.GLCM_Homogeneity1 = mean(Radiomics.GLCM_Homogeneity1);
Radiomics.GLCM_Homogeneity2 = mean(Radiomics.GLCM_Homogeneity2);
Radiomics.GLCM_IMC1 = mean(Radiomics.GLCM_IMC1);
Radiomics.GLCM_IMC2 = mean(Radiomics.GLCM_IMC2);
Radiomics.GLCM_IDN = mean(Radiomics.GLCM_IDN);
Radiomics.GLCM_IDMN = mean(Radiomics.GLCM_IDMN);
Radiomics.GLCM_InverseVariance = mean(Radiomics.GLCM_InverseVariance);
Radiomics.GLCM_MaximumProbability = mean(Radiomics.GLCM_MaximumProbability);
Radiomics.GLCM_SumAverage = mean(Radiomics.GLCM_SumAverage);
Radiomics.GLCM_SumEntropy = mean(Radiomics.GLCM_SumEntropy);
Radiomics.GLCM_SumVariance = mean(Radiomics.GLCM_SumVariance);
Radiomics.GLCM_Variance = mean(Radiomics.GLCM_Variance);