function excludedTrials = rd_temporalAttentionExcludeTrials(expt, ...
    startingTrial, nGoodTrialsToExclude)
%
% function excludeTrials = rd_temporalAttentionExcludeTrials(expt, startingTrial, nGoodTrialsToExclude)

if nargin==1
    startingTrial = 1;
    nGoodTrialsToExclude = 64;
end

eye = expt.eye;
eyeGood = eye.fixCue & eye.fixT1 & eye.fixT2;

excludeTrials = zeros(size(eyeGood));
nTrials = numel(excludeTrials);

excludedTrialsCounter = 0;
for iTrial = startingTrial:nTrials
    if eyeGood(iTrial)
        excludeTrials(iTrial) = 1;
        excludedTrialsCounter = excludedTrialsCounter + 1;
    end
    if excludedTrialsCounter==nGoodTrialsToExclude
        break
    end
end

excludeTrials = logical(excludeTrials);
excludedTrials = expt.trialOrder(excludeTrials);

