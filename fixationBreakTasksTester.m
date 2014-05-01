% fixationBreakTasksTester.m

% dummy vars
window = 0;
p = 0;

trials = fullfact([4 4]);
nTrials = size(trials,1);

trialOrder = randperm(nTrials)';
fixations = mod(trialOrder,2);

iTrial = 1;
while iTrial <= nTrials
    iTrial
    trialIdx = trialOrder(iTrial);
    
    if iTrial <= length(fixations)
        fixation = fixations(iTrial);
    else
        fixation = 1;
    end
    
    [stopThisTrial trialOrder, nTrials] = fixationBreakTasks(...
        fixation, window, p, trialOrder, iTrial, nTrials);
    iTrial = iTrial+1;
    
    if stopThisTrial
        continue
    end
    
    trials(trialIdx,3) = fixation;
end