% fixADPilotRun1.m

%% Load the saved workspace
load([pathToExpt('data') '/adPilot_cb_tc64-100_soa1000-1250_run01_WORKSPACE.mat'])

%% Set missing data to nan
skippedTrials = trials(:,10)==0;
trials(skippedTrials, 6:end) = NaN;

%% Store expt info
expt.subjectID = subjectID;
expt.p = p;
expt.timing = timing;
expt.trialOrder = trialOrder;
expt.trials_headers = trials_headers;
expt.trials = trials;
expt.targetRotations = targetRotations;
expt.targetPhases = targetPhases;
expt.trialsPresented = trialsPresented;

%% Calculate more timing things
% note these values will sometimes be wrong, and we aren't bothering to
% correct them
expt.timing.dur.im1 = expt.timing.timeBlank1 - expt.timing.timeIm1;
expt.timing.dur.im2 = expt.timing.timeBlank2 - expt.timing.timeIm2;
expt.timing.dur.cueIm1SOA = expt.timing.timeIm1 - expt.timing.timeCue;
expt.timing.dur.cueIm2SOA = expt.timing.timeIm2 - expt.timing.timeCue;
expt.timing.dur.im1Im2SOA = expt.timing.timeIm2 - expt.timing.timeIm1;

%% Analyze and save data
% [expt results] = rd_analyzeTemporalAttention(expt, saveData, saveFigs, plotTimingFigs, saveTimingFigs);
[expt results] = rd_analyzeTemporalAttention(expt, 1, 1, 0, 0);