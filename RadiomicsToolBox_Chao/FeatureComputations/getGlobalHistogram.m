function [HistFeas] = getGlobalHistogram(ROIonly, N)

% PRELIMINARY
Nbins=N;
vectorValid = ROIonly(~isnan(ROIonly));
numVoxels = numel(vectorValid);
histo = hist(vectorValid,Nbins);
histoPercent = histo./(sum(histo(:)));
vectNg = 1:Nbins;
u = histoPercent*vectNg';


% COMPUTATION OF HistFeas
% 1. Variance
variance = 0;
for i=1:Nbins
    variance = variance+histoPercent(i)*(i-u)^2;
end
sigma = sqrt(variance);
HistFeas.globalHistogram_Variance = variance;

% 2. Skewness
skewness = 0;
for i = 1:Nbins
    skewness = skewness+histoPercent(i)*(i-u)^3;
end
skewness = skewness/sigma^3;
HistFeas.globalHistogram_Skewness = skewness;

% 3. Kurtosis
kurtosis = 0;
for i = 1:Nbins
    kurtosis = kurtosis+histoPercent(i)*(i-u)^4;
end
kurtosis = (kurtosis/sigma^4) - 3;
HistFeas.globalHistogram_Kurtosis = kurtosis;

% 4. Entropy
%HistFeas.globalHistogram_Entropy=-sum(histoPercent.*log2(histoPercent+realmin));
HistFeas.globalHistogram_Entropy=entropy(histoPercent);

% 5. Uniformity
HistFeas.globalHistogram_Uniformity=sum(histoPercent.*histoPercent);

% 6 Mean 
HistFeas.globalHistogram_Mean = mean(histoPercent);

% 7 Median
HistFeas.globalHistogram_Median = median(histoPercent);

%8 Range
HistFeas.globalHistogram_Range = range(histoPercent);%difference between large and min of histogram.

%9 Interquartile range
p75 = prctile(histoPercent,75);
p25 = prctile(histoPercent,25);
HistFeas.globalHistogram_InterquartileRange = p75 - p25;

%10 quartile coefficient of dispersion
HistFeas.globalHistogram_QuartileCoefficientDispersion = (p75-p25)/(p75+p25);

% 11
HistFeas.globalHistogram_Percentile10 = prctile(histoPercent,10);

% 12
HistFeas.globalHistogram_Percentile90 = prctile(histoPercent,90);

% 13
HistFeas.globalHistogram_Mode = mode(histoPercent);

% 14
HistFeas.globalHistogram_MeanAbsoluteDeviation = mad(histoPercent);

%15 robust mean absolute deviation
p10 = prctile(histoPercent,10);
p90 = prctile(histoPercent,90);
histoPercent10_90 = histoPercent(histoPercent >= p10 & histoPercent <= p90);
HistFeas.globalHistogram_RobustMeanAbsoluteDeviation = mad(histoPercent10_90);

%16 median absolute deviation
HistFeas.globalHistogram_MedianAbsoluteDeviation = sum(abs(histoPercent-median(histoPercent)))/numel(histoPercent);

%17 coefficient of variation
HistFeas.globalHistogram_CoefficientVariation = sigma/u;

end