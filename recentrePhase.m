function phaseData = recentrePhase(phaseData, phaseCentre)
% phaseData = recentrePhase(phaseData, phaseCentre)
%
% Function converts phase values to values within a chosen phase range
% centred around input variable phaseCentre: [phaseCentre-pi
% phaseCentre+pi].
% Input: phaseData - a vector or a matrix of phase values.
%        phaseCentre - the centre of the new phase range.
% Output: phaseData - recentred phase data.

upperLimit = phaseCentre + pi;
lowerLimit = phaseCentre - pi;
for i = 1:numel(phaseData)
  if phaseData(i) > upperLimit
    while phaseData(i) > upperLimit
      phaseData(i) = phaseData(i) - 2*pi;
    end
  elseif phaseData(i) < lowerLimit
    while phaseData(i) < lowerLimit
      phaseData(i) = phaseData(i) + 2*pi;
    end
  end
end