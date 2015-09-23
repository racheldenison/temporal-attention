function tcSame = rd_findOneBackSameDiff(expt)

% load([pathToExpt('data') '/E4_contrast_cbD10/ax/ax_cbD10_tilt*_tc16-64_soa1000-1300_run08_TemporalAttention_T1T2all_20150923.mat'])


%% setup
trialOrder = expt.trialOrder;
trialOrderRun = expt.trialOrderRun;
fixation = expt.fixation;
trials_headers = expt.trials_headers;
trials = expt.trials;

targetContrastIdx = strcmp(trials_headers,'targetContrast');

targetContrasts = trials(:,targetContrastIdx);
nTrials = size(trials,1);

%% find whether target contrast was the same or different as previous trial
runs = unique(trialOrderRun);
nRuns = numel(runs);
for iRun = 1:nRuns
    tr = (1:nTrials/nRuns) + (nTrials/nRuns)*(iRun-1);
    tc = targetContrasts(tr);
    to = trialOrder(trialOrderRun==runs(iRun));
    f = fixation(trialOrderRun==runs(iRun));
    tcso = nan(size(to)); % in order of presentation
    tcs = nan(size(tr)); % order of trials matrix
    
    for iTrial = 2:numel(to)
        % don't look at trials where there was a fixation break or
        % trials where the previous trial had a fixation break
        if f(iTrial)==1 && f(iTrial-1)==1
            trialIdx = to(iTrial); % index into trials matrix
            trialIdx1 = to(iTrial-1);
            tcso(iTrial) = tc(trialIdx)==tc(trialIdx1);
            tcs(trialIdx) = tc(trialIdx)==tc(trialIdx1);
        end
    end
    
    tcSame(tr) = tcs;
end
   
tcSame = tcSame';


