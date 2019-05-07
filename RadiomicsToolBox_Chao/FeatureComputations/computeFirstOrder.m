function [Radiomics] = computeFirstOrder(I)
% I is ROIonly

m = floor(min(min(min(I))));
M = ceil(max(max(max(I))));
levels = m:M;

X = I(~isnan(I));
N = numel(X);
% X = reshape(I, [N 1]); % deleted cuz useless with "I = I(~isnan(I))"
X_Bar = sum(X)/N;
if max(levels) == 0 && min(levels) == 0
    [P, ~] = histcounts(X);
else
    [P, ~] = histcounts(X, levels);
end
P = P/N;
Radiomics.firstOrder_Energy = sum(P.^2);
Radiomics.firstOrder_Entropy = -P(P~=0)*log2(P(P~=0)');
Radiomics.firstOrder_Kurtosis = kurtosis(X, 1);% Kurtosis is a measure of how outlier-prone a distribution is. The kurtosis
% of the normal distribution is 3. Distributions that are more outlier-prone
% than the normal distribution have kurtosis greater than 3; distributions
% that are less outlier-prone have kurtosis less than 3.

Radiomics.firstOrder_Maximum = max(X);
Radiomics.firstOrder_Mean = X_Bar;


%   Mean deviation (also called mean absolute deviation)
Radiomics.firstOrder_MeanAbsoluteDeviation = mad(X);

Radiomics.firstOrder_Median = median(X);
Radiomics.firstOrder_Minimum = min(X);
Radiomics.firstOrder_Range = range(X);


%   Root mean square (RMS)
Radiomics.firstOrder_RootMeanSquare = rms(X);

% Skewness is a measure of the asymmetry of the data around the sample mean.
% If skewness is negative, the data are spread out more to the left of the mean
% than to the right. If skewness is positive, the data are spread out more to the
% right. The skewness of the normal distribution (or any perfectly symmetric
% distribution) is zero.
Radiomics.firstOrder_Skewness = skewness(X, 1);

Radiomics.firstOrder_StandardDeviation = std(X);
Radiomics.firstOrder_Uniformity = P*P';
Radiomics.firstOrder_Variance = var(X);

% features above are from "Chihyanghsu0805/Radiomics".

% features below are from "adityaapte/CERR".

% Median absolute deviation
Radiomics.firstOrder_MedianAbsoluteDeviation = sum(abs(X-Radiomics.firstOrder_Median))/numel(X);

%   P10
p10 = prctile(X,10);
Radiomics.firstOrder_Percentile10 = p10;

%   P90
p90 = prctile(X,90);
Radiomics.firstOrder_Percentile90 = p90;

X10_90 = X(X >= p10 & X <= p90);

%   Robust Mean Absolute Deviation
Radiomics.firstOrder_RobustMeanAbsoluteDeviation  = mad(X10_90);

% Inter-Quartile Range (IQR)
% P75 - P25
p75 = prctile(X,75);
p25 = prctile(X,25);
Radiomics.firstOrder_InterquartileRange = p75 - p25;

% Quartile coefficient of Dispersion
Radiomics.firstOrder_QuartileCoefficientDispersion = (p75-p25)/(p75+p25);

% Coefficient of variation
Radiomics.firstOrder_CoefficientVariation = Radiomics.firstOrder_StandardDeviation / Radiomics.firstOrder_Mean;
