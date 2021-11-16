function [U2, pval] = WatsonU2_goodnessOfFit(F, sumTerms)
% [U2, pval] = WatsonU2_goodnessOfFit(F, sumTerms)
%
% Function computes Watson U2 goodness-of-fir statistic.
% Input: F - fitted cummulative distribution fucntion of your fit.
%        sumTerms - number of exponential sum terms in the estimation of
%                   p-value. The default value is 30 as a larger number of
%                   terms provides negligible improvement in convergence.
% Output: U2 - Watson U2 goodness-of-fit statistic.
%         pval - associated p-value. Larger U2 values give smaller
%                p-values. Therefore, the bigger is the U2 statistic, the
%                poorer the fit. If your aim is to test the similarity of
%                the fitted and empirical distributions, you should,
%                therefore, invert your p-value:
%                new p-value = 1 - original p-value.

% by Martynas Dervinis (martynas.dervinis@gmail.com)

if nargin < 2
  sumTerms = 30;
end

n = numel(F);
u_bar = mean(F);

% Rank the sample
[~, sortOrder] = sort(F);
ranks = zeros(1,n);
for iF = 1:n
  ranks(sortOrder(iF)) = iF;
end

% Calculate U2 goodness-of-fit statistic
U2 = 0;
for iF = 1:n
  U2 = U2 + (F(iF) - (u_bar - 0.5) - ((2*ranks(iF) - 1)/(2*n)))^2;
end
U2 = U2 + (1/(12*n));

% Estimate the p-value
pval = 0;
for iSum = 1:sumTerms
  pval = pval + ((-1)^(iSum-1))*2*exp(-2*(iSum^2)*(pi^2)*U2);
end