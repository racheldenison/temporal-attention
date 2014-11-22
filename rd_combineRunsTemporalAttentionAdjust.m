% rd_combineRunsTemporalAttentionAdjust.m

%% setup
subject = 'bl_a1_tc100_soa1000-1250';
runs = 1:4;
combinedRun = 9;
nRuns = numel(runs);

saveData = 1;
saveFigs = 1;

expName = 'E3_adjust';
% dataDir = 'data';
% figDir = 'figures';
dataDir = pathToExpt('data');
figDir = pathToExpt('figures');
dataDir = sprintf('%s/%s/%s', dataDir, expName, subject(1:2));
figDir = sprintf('%s/%s/%s', figDir, expName, subject(1:2));

%% initializations
subjectID = sprintf('%s_run%02d', subject, combinedRun);
trials = [];
targetRotations = [];
trialOrder = [];

%% get data from each run
for iRun = 1:nRuns
    run = runs(iRun);
    dataFile = dir(sprintf('%s/%s_run%02d*', dataDir, subject, run));
    data = load(sprintf('%s/%s', dataDir, dataFile.name));
    
    trials = [trials; data.expt.trials];
    targetRotations = [targetRotations; data.expt.targetRotations];
    trialOrder = [trialOrder; data.expt.trialOrder + size(trialOrder,2)];
    
    if strcmp(subjectID(1:3),'xx*')
        timing = [];
    else
        timing(iRun) = data.expt.timing;
    end
end

%% make the combined expt
expt.subjectID = subjectID;
expt.runs = runs;
expt.p = data.expt.p; % assume all runs identical
expt.timing = timing;
expt.trialOrder = trialOrder;
expt.trials_headers = data.expt.trials_headers; % assume all runs identical
expt.trials = trials;
expt.targetRotations = targetRotations;

%% analyze data
[expt results] = rd_analyzeTemporalAttentionAdjust(expt);

%% save data
% saving data and figs separately in order to save them into the mcq
% directory, and not locally
if saveData
    fileName = sprintf('%s/%s_TemporalAttentionAdjust_T1T2all_%s.mat', dataDir, subjectID, datestr(now, 'yyyymmdd'));
    save(fileName, 'expt', 'results')
end

%% save figs
if saveFigs
    figNames = {'error','acc','rt'};
    rd_saveAllFigs([], figNames, [subjectID '_TemporalAttentionAdjust_T1T2all'], figDir)
end
    