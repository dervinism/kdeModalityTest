function [resultsTable, modalityFits, testStatistic, pValDiff, pValFit, modes, fits, fitPts] = kdeModalityTest(alpha, type, lambda_0, m_0, modeStop, binSize)
% [resultsTable, modalityFits, testStatistic, pValDiff, pValFit, modes, fits, fitPts] = kdeModalityTest(alpha, type, lambda_0, m_0, modeStop, binSize)
%
% Function determines data modality (number of modes in the data) and data
% mode properties. Modality estimation is based on kernel density excess
% mass estimate approach described by Fisher and Marron (2001). The linear
% goodness-of-fit test is based on Cramer-von Misses T_k statistic (not
% implemented yet). The circular goodness-of-fit test is based on Watson U2
% statistic as in Watson (1961).
%
% Input: alpha - data sample vector. If analysis type is circular, data
%                should be in radians.
%        type - data type vector: 'linear' or 'circular'. Default is
%               'linear'.
%        lambda_0 - modality estimation parameter for weeding out low
%                   amplitude modes. Default is 0.
%        m_0 - modality estimation parameter for weeding out low mass
%              modes. m_0 must be within the range of 0.01-1. Default is
%              0.05 (corresponds to 5% of the overall distribution).
%        modeStop - the maximum number of modes to consider. The default
%                   value is inf but in reality it is limited by the
%                   fitting sampling size which is 1/500th of the data
%                   range if data is linear and 1/(2*pi) if data is
%                   circular.
%        binSize - a two-entry vector describing the bin size. First entry
%                  describes the bin size for mode estimation part. The
%                  second entry describes the bin size for statistical
%                  testing part. Default size for the first entry is
%                  1/500th for linear data and 1/(2*pi) for circular data.
%                  Default size for the second entry is 1/10th of the
%                  first entry.
%
% Output: resultsTable is the main output variable that will give you a
%                     comprehensive summary of data modality. It is a
%                     matrix with rows corresponding to different data
%                     descriptors. Row one indicates the number of modes.
%                     Row two is the bandwidth. Row three is the excess
%                     mass estimate. Row four is the goodness-of-fit
%                     statistic (linear T_k or circular U2). Row five is a
%                     corresponding p-value for the difference between the
%                     fitted modal distribution and the actual data. If
%                     this p-value exceeds 0.05, that means that the null
%                     hypothesis of the fitted modal and empirical data
%                     distributions being drawn from the same distribution
%                     can no longer be rejected. The row six is the
%                     corresponding p-value for the similarity of the
%                     fitted and empirical distributions:
%                     p-value_similarity = 1 - p-value_difference. Fits
%                     should have similarity p-value of less than 0.05 for
%                     them to be deemed as good enough. You should use the
%                     cut-off of 0.05 for deciding the minimal modality
%                     (lower limit) of your data. Increasing the number of
%                     modes in your fitted distribution is most likely to
%                     provide better fits but these are also likely to be
%                     overfits.
%         modalityFits - rows one to three of resultsTable.
%         testStatistic - row four of resultsTable.
%         pValDiff - row five of resultsTable.
%         pValFit - row six of resultsTable.
%         modes - a cell array with fields holding mode properties
%                 corresponding to columns of resultsTable. Fields have the
%                 following properties:
%                   peaks - mode peak locations in linear units or radians.
%                   centreMass - mode mass centres.
%                   masses - mode masses in terms of cumulitive
%                            distribution (<= 1). This property could be
%                            thought of as mode size in terms of proportion
%                            of the total sample.
%         fits - a matrix with rows corresponding to function fits of
%                various modalities (columns of resultsTable).
%         fitPts - fit points. Use plot(fitPts,fits(1-to-n,:)) to display
%                  the fitted modal distribution function.
%
% Dependencies: circ_ksdensity (Muir, 2020) and WatsonU2_goodnessOfFit,
%               recentrePhase.
%
% References: N. I. Fisher and J. S. Marron, "Mode Testing via the Excess
%               Mass Estimate," Biometrika 88 (2), 499-517 (2001).
%             G. S. Watson, “Goodness-of-Fit Tests on a Circle. I,”
%               Biometrika 48 (1–2), 109–114 (1961).
%             D. Muir, "Kernel density estimation for circular functions"
%               (https://www.mathworks.com/matlabcentral/fileexchange/44072
%               -kernel-density-estimation-for-circular-functions), MATLAB
%               Central File Exchange. Retrieved June 7, 2020.

% by Martynas Dervinis (martynas.dervinis@gmail.com)

% Parameters of the Excess Mass Estimate algorithm
if nargin < 5 || isempty(modeStop)
  modeStop = inf;
end
if nargin < 4 || isempty(m_0)
  m_0 = 0.05;
elseif m_0 < 0.01
  m_0 = 0.01;
elseif m_0 > 1
  m_0 = 1;
end
if nargin < 3 || isempty(lambda_0)
  lambda_0 = 0;
end
if nargin < 2 || isempty(type) || strcmpi(type, '')
  type = 'linear';
elseif ~strcmpi(type, 'linear') && ~strcmpi(type, 'circular')
  error('The function only works with linear or circular data types');
end
if strcmpi(type, 'linear')
  domain = [min(alpha) max(alpha)];
elseif strcmpi(type, 'circular')
  domain = [-pi pi];
end
if nargin < 6 || isempty(binSize)
  if strcmpi(type, 'linear')
    nSamples = 500;
    binSize(1) = (domain(2) - domain(1))/nSamples;
  elseif strcmpi(type, 'circular')
    nSamples = 360;
    binSize(1) = (2*pi)/nSamples;
  end
  binSize(2) = binSize(1)/10;
end
samplePts = domain(1)+binSize(1):binSize(1):domain(2);
denseSamplePts = domain(1)+binSize(2):binSize(2):domain(2);
if ~(nargin < 6 || isempty(binSize))
  nSamples = numel(samplePts);
end

% Preprocess data
alpha = alpha(~isnan(alpha));
if strcmpi(type, 'circular')
  alpha = recentrePhase(alpha, 0);
end

% Estimate the minimum number of modes required to describe the data
resultsTable = []; modalityFits = []; testStatistic = []; pValDiff = []; pValFit = []; modes = []; fits = []; fitPts = [];
h = ((nSamples/4)*binSize(1)):-binSize(1):binSize(1);
for iH = 1:numel(h) % loop through increasingly smaller h size
%   disp([num2str(iH) '/' num2str(numel(h))]);
  
  % Step 1: Kernel density estimation
  if strcmpi(type, 'linear')
    f = ksdensity(alpha, samplePts, 'Bandwidth',h(iH));
  elseif strcmpi(type, 'circular')
    f = circ_ksdensity(alpha, samplePts, domain, h(iH)); % The kernel function is wrapped Gaussian, not von Mises
    extendedSamplePts = [samplePts-2*pi samplePts samplePts+2*pi]; %#ok<*NASGU>
  end
  cdf = cumsum(f)./sum(f);
  
  % Step 2: Maxima localisation and their base estimation
  %figure; plot(samplePts, f); hold on; plot(samplePts, 1./f); hold off;
  if strcmpi(type, 'linear')
    extendedF = f;
    extendedCDF = cdf;
  elseif strcmpi(type, 'circular')
    extendedF = [f f f];
    extendedCDF = [cdf 1+cdf 2+cdf];
  end
  [peaks, peakLocs] = findpeaks(extendedF, 'Threshold',lambda_0);
  [~, troughLocs] = findpeaks(1./extendedF);
  troughs = extendedF(troughLocs);
  if peakLocs(1) ~= 1 && peakLocs(1) < troughLocs(1)
    troughLocs = [1 troughLocs];
    troughs = [extendedF(1) troughs];
  end
  if peakLocs(end) ~= numel(extendedF) && peakLocs(end) > troughLocs(end)
    troughLocs = [troughLocs numel(extendedF)];
    troughs = [troughs extendedF(end)];
  end
  if troughLocs(1) ~= 1 && troughLocs(1) < peakLocs(1)
    peakLocs = [1 peakLocs];
    peaks = [extendedF(1) peaks];
  end
  if troughLocs(end) ~= numel(extendedF) && troughLocs(end) > peakLocs(end)
    peakLocs = [peakLocs numel(extendedF)];
    peaks = [peaks extendedF(end)];
  end
  E = zeros(1,numel(peaks));
  centreMass = zeros(1,numel(peaks));
  base = {};
  baseLocs = {};
  lambda = {};
  lambdaLocs = {};
  for iPks = 1:numel(peaks)
    if iPks == 1 && peakLocs(iPks) == 1
      base{iPks} = [troughs(1) troughs(1)]; %#ok<*AGROW>
      baseLocs{iPks} = [1 troughLocs(1)];
    elseif iPks == numel(peaks) && peakLocs(iPks) == numel(extendedF)
      base{iPks} = [troughs(end) troughs(end)];
      baseLocs{iPks} = [troughLocs(end) numel(extendedF)];
    else
      [~, troughOnRight] = find(troughLocs > peakLocs(iPks), 1);
      baseLocs{iPks} = [troughLocs(troughOnRight-1) troughLocs(troughOnRight)];
      base{iPks} = [troughs(troughOnRight-1) troughs(troughOnRight)];
    end
    [lambda{iPks}, lambdaLocs{iPks}] = baseEstimation(extendedF(baseLocs{iPks}(1):baseLocs{iPks}(2)), base{iPks}, baseLocs{iPks}, lambda_0);
    
    % Step 3: Calculation of excess mass associated with each maximum
    [E(iPks), centreMass(iPks)] = excessMass(extendedF, extendedCDF, lambda{iPks}, lambdaLocs{iPks});
  end
  
%   figure; hold on
%   plot(extendedSamplePts,extendedF);
%   plot(extendedSamplePts(peakLocs),peaks, 'g.', 'MarkerSize',10);
%   plot(extendedSamplePts(troughLocs),troughs, 'r.', 'MarkerSize',10);
%   hold off
  
  % Step 4: Elimination of isolated minor modes
  minorModes = logical(E <= m_0);
  if sum(minorModes)
    modes2eliminate = [];
    for iMode = 1:numel(E)
      if minorModes(iMode) && max(base{iMode}) <= lambda_0
        modes2eliminate = [modes2eliminate iMode];
      end
    end
    if ~isempty(modes2eliminate)
      E = E(~modes2eliminate);
      centreMass = centreMass(~modes2eliminate);
      peaks = peaks(~modes2eliminate);
      peakLocs = peakLocs(~modes2eliminate);
      base = base(~modes2eliminate);
      baseLocs = baseLocs(~modes2eliminate);
      lambda = lambda(~modes2eliminate);
      lambdaLocs = lambdaLocs(~modes2eliminate);
    end
  end
  
%   minorModes = logical(E <= m_0);
%   hold on
%   plot(extendedSamplePts(peakLocs(minorModes)),peaks(minorModes), 'rx', 'MarkerSize',10)
%   hold off
  
  % Step 5: Merger of the remaining minor modes
  iMode = 0;
  while true
    minorModes = logical(E <= m_0);
    if sum(minorModes) && iMode == numel(E) && ~eliminate
      break
    elseif sum(minorModes)
      for iMode = 1:numel(E)
        if minorModes(iMode)
          [~, ind] = max(base{iMode});
          eliminate = false;
          if iMode == 1 && numel(E) > 1 && min(base{iMode}) > min(base{iMode+1})
            peaks2merge = [peaks(iMode) peaks(iMode+1)];
            peakLocs2merge = [peakLocs(iMode) peakLocs(iMode+1)];
            [peaks(iMode+1), whichPeak] = max(peaks2merge);
            peakLocs(iMode+1) = peakLocs2merge(whichPeak);
            base{iMode+1} = [base{iMode}(1) base{iMode+1}(2)];
            baseLocs{iMode+1} = [baseLocs{iMode}(1) baseLocs{iMode+1}(2)];
            [lambda{iMode+1}, lambdaLocs{iMode+1}] = baseEstimation(extendedF(baseLocs{iMode+1}(1):baseLocs{iMode+1}(2)),...
              base{iMode+1}, baseLocs{iMode+1}, lambda_0);
            [E(iMode+1), centreMass(iMode+1)] = excessMass(extendedF, extendedCDF, lambda{iMode+1}, lambdaLocs{iMode+1});
            eliminate = true;
          elseif (iMode == numel(E) && numel(E) > 1 && min(base{iMode}) > min(base{iMode-1})) || (ind == 1 && iMode ~= 1)
            peaks2merge = [peaks(iMode-1) peaks(iMode)];
            peakLocs2merge = [peakLocs(iMode-1) peakLocs(iMode)];
            [peaks(iMode-1), whichPeak] = max(peaks2merge);
            peakLocs(iMode-1) = peakLocs2merge(whichPeak);
            base{iMode-1} = [base{iMode-1}(1) base{iMode}(2)];
            baseLocs{iMode-1} = [baseLocs{iMode-1}(1) baseLocs{iMode}(2)];
            [lambda{iMode-1}, lambdaLocs{iMode-1}] = baseEstimation(extendedF(baseLocs{iMode-1}(1):baseLocs{iMode-1}(2)),...
              base{iMode-1}, baseLocs{iMode-1}, lambda_0);
            [E(iMode-1), centreMass(iMode-1)] = excessMass(extendedF, extendedCDF, lambda{iMode-1}, lambdaLocs{iMode-1});
            eliminate = true;
          elseif (ind == 2 && iMode ~= numel(E))
            peaks2merge = [peaks(iMode) peaks(iMode+1)];
            peakLocs2merge = [peakLocs(iMode) peakLocs(iMode+1)];
            [peaks(iMode+1), whichPeak] = max(peaks2merge);
            peakLocs(iMode+1) = peakLocs2merge(whichPeak);
            base{iMode+1} = [base{iMode}(1) base{iMode+1}(2)];
            baseLocs{iMode+1} = [baseLocs{iMode}(1) baseLocs{iMode+1}(2)];
            [lambda{iMode+1}, lambdaLocs{iMode+1}] = baseEstimation(extendedF(baseLocs{iMode+1}(1):baseLocs{iMode+1}(2)),...
              base{iMode+1}, baseLocs{iMode+1}, lambda_0);
            [E(iMode+1), centreMass(iMode+1)] = excessMass(extendedF, extendedCDF, lambda{iMode+1}, lambdaLocs{iMode+1});
            eliminate = true;
          end
          if eliminate
            peaks(iMode) = [];
            peakLocs(iMode) = [];
            trough2eliminate = logical(troughLocs > peakLocs2merge(1) & troughLocs < peakLocs2merge(2));
            troughs(trough2eliminate) = [];
            troughLocs(trough2eliminate) = [];
            base(iMode) = [];
            baseLocs(iMode) = [];
            lambda(iMode) = [];
            lambdaLocs(iMode) = [];
            E(iMode) = [];
            centreMass(iMode) = [];
            break
          end
        end
      end
    else
      if ~isempty(E)
        if (isempty(troughLocs) || peakLocs(1) < troughLocs(1)) && peakLocs(1) ~= 1
          [firstTrough, firstTroughLoc] = min(extendedF(1:peakLocs(1)));
          troughs = [firstTrough troughs];
          troughLocs = [firstTroughLoc troughLocs];
        elseif (isempty(troughLocs) || peakLocs(end) > troughLocs(end)) && peakLocs(end) ~= numel(extendedF)
          [lastTrough, lastTroughLoc] = min(extendedF(peakLocs(end):end));
          troughs = [troughs lastTrough];
          troughLocs = [troughLocs lastTroughLoc];
        end
      end
      break
    end
  end
  
  % Step 6: Elimination of minor modes that couldn't be merged (existing at data sample ends)
  minorModes = logical(E <= m_0);
  if sum(minorModes)
    E(minorModes) = [];
    centreMass(minorModes) = [];
    peaks(minorModes) = [];
    peakLocs(minorModes) = [];
    base(minorModes) = [];
    baseLocs(minorModes) = [];
    lambda(minorModes) = [];
    lambdaLocs(minorModes) = [];
  end
  
%   figure; hold on
%   plot(extendedSamplePts,extendedF);
%   plot(extendedSamplePts(peakLocs),peaks, 'g.', 'MarkerSize',10);
%   plot(extendedSamplePts(troughLocs),troughs, 'r.', 'MarkerSize',10);
%   hold off
  
  % Step 7: Estimation of the excess mass
  if strcmpi(type, 'circular')
    E = E(peakLocs > numel(f) & peakLocs <= 2*numel(f));
    centreMass = centreMass(peakLocs > numel(f) & peakLocs <= 2*numel(f));
    peaks = peaks(peakLocs > numel(f) & peakLocs <= 2*numel(f));
    peakLocs = peakLocs(peakLocs > numel(f) & peakLocs <= 2*numel(f)) - numel(f);
    troughs = troughs(troughLocs > numel(f) & troughLocs <= 2*numel(f));
    troughLocs = troughLocs(troughLocs > numel(f) & troughLocs <= 2*numel(f)) - numel(f);
  end
  sortedE = sort(E, 'descend');
  if numel(sortedE) > 1
    S_k = sum(sortedE(2:end));
  else
    S_k = 0;
  end
  
  if numel(E) > modeStop
    return % This is where the code stops if the maximum number of modes has been exceeded
  end
  
%   figH = figure; hold on
%   plot(samplePts,f);
%   plot(samplePts(peakLocs),peaks, 'g.', 'MarkerSize',10);
%   plot(samplePts(troughLocs),troughs, 'r.', 'MarkerSize',10);
%   plot(extendedSamplePts(centreMass),extendedF(centreMass), 'b.', 'MarkerSize',10);
%   plot(samplePts,circ_ksdensity(alpha, samplePts, domain, h(end)));
%   hold off
%   disp([numel(E), S_k]);
%   close all
  
  % Step 8: Calculation of appropriate statistics
  if ~isempty(E)
    if strcmpi(type, 'linear') % Currently not implemented
      F_alpha = ksdensity(alpha, alpha, 'Bandwidth',h(iH), 'Function','cdf');
      %stat = CramervonMises(F_alpha);
      stat = NaN;
      pval = NaN;
    elseif strcmpi(type, 'circular')
      f_alpha = circ_ksdensity(alpha, denseSamplePts, domain, h(iH));
      F = cumsum(f_alpha)./sum(f_alpha);
      F_alpha = zeros(1,numel(alpha));
      for iF = 1:numel(alpha)
        [~, closestInd] = min(abs(denseSamplePts - alpha(iF)));
        F_alpha(iF) = F(closestInd);
      end
      [stat, pval] = WatsonU2_goodnessOfFit(F_alpha);
%       disp([stat, pval]);
%       close(figH);
    end
  else
    stat = NaN;
    pval = NaN;
  end
  
  % Step 9: Compilation of output variables
  if ~isempty(E)
    if isempty(modalityFits)
      modalityFits = [numel(E); h(iH); S_k];
      testStatistic = stat;
      pValDiff = min([pval 1]);
      pValFit = max([1 - pval 0]);
      resultsTable = [numel(E); h(iH); S_k; stat; pValDiff; pValFit];
      modes{1}.peaks = samplePts(peakLocs);
      modes{1}.centreMass = recentrePhase(extendedSamplePts(centreMass), 0);
      modes{1}.masses = E;
      fits = f;
      fitPts = samplePts;
    elseif (resultsTable(1,end) == numel(E) && resultsTable(5,end) < pval) ||...
        (resultsTable(1,end) == numel(E) && resultsTable(5,end) == pval && resultsTable(2,end) < S_k)
      modalityFits(:,end) = [numel(E); numel(E); S_k];
      testStatistic(end) = stat;
      pValDiff(end) = min([pval 1]);
      pValFit(end) = max([1 - pval 0]);
      resultsTable(:,end) = [numel(E); h(iH); S_k; stat; min([pval 1]); max([1 - pval 0])];
      modes{end}.peaks = samplePts(peakLocs);
      modes{end}.centreMass = recentrePhase(extendedSamplePts(centreMass), 0);
      modes{end}.masses = E;
      fits(end,:) = f;
    elseif resultsTable(1,end) < numel(E)
      modalityFits = [modalityFits [numel(E); h(iH); S_k]];
      testStatistic = [testStatistic stat];
      pValDiff = [pValDiff min([pval 1])];
      pValFit = [pValFit max([1 - pval 0])];
      resultsTable = [resultsTable [numel(E); h(iH); S_k; stat; min([pval 1]); max([1 - pval 0])]];
      modes{numel(modes)+1}.peaks = samplePts(peakLocs);
      modes{end}.centreMass = recentrePhase(extendedSamplePts(centreMass), 0);
      modes{end}.masses = E;
      fits = [fits; f];
    end
  end
end


%% Local functions
function [baseEst, baseLocsEst] = baseEstimation(f, baseInit, baseInitLocs, lambda_0)
% Function estimates mode's base (lambda) and base locations.

if baseInitLocs(1) == 1 || baseInitLocs(2) == numel(f)
  thresholdedBase = max([min(baseInit) lambda_0]);
else
  thresholdedBase = max([baseInit lambda_0]);
end
peakWithBaseInds = logical(f >= thresholdedBase);
peakWithBase = f(peakWithBaseInds);
baseEst = [min(peakWithBase) min(peakWithBase)];

peakWithBaseLocsInit = baseInitLocs(1):baseInitLocs(2);
peakWithBaseLocs = peakWithBaseLocsInit(peakWithBaseInds);
baseLocsEst = [peakWithBaseLocs(1) peakWithBaseLocs(end)];


function [E, centreMass] = excessMass(pdf, cdf, lambda, lambdaLocs)
% Function extimates mode's mass and the centre of mass

% Excess mass
AUC = sum(pdf(lambdaLocs(1):lambdaLocs(2)));
ABC = AUC - lambda(1)*numel(pdf(lambdaLocs(1):lambdaLocs(2)));
if lambdaLocs(1) == 1
  cdfOI = cdf(lambdaLocs(2));
  E = cdfOI*(ABC/AUC);
else
  cdfOI = cdf(lambdaLocs(2)) - cdf(lambdaLocs(1)-1);
  E = cdfOI*(ABC/AUC);
end

% Mass centre
cdfRangeOI = cdf(lambdaLocs(1):lambdaLocs(2));
if lambdaLocs(1) == 1
  centreMass = cdfOI/2;
else
  centreMass = cdf(lambdaLocs(1)-1) + cdfOI/2;
end
[~, centreMass] = min(abs(cdfRangeOI - centreMass));
centreMass = lambdaLocs(1) + centreMass - 1;