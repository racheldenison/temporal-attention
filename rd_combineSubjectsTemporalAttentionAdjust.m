% rd_combineSubjectsTemporalAttentionAdjust.m

%% setup
subjectIDs = {'bl','rd','id','ec','ld','en','sj','ml','ca','jl','ew','jx'};
nSubjects = numel(subjectIDs);
run = 9;
expName = 'E3_adjust';

%% initialize super-subject structure
SS.expt.trials = [];
SS.expt.targetRotations = [];
for iV = 1:3
    for iT = 1:2
        SS.results.totals.all{iV,iT} = [];
    end
end

%% combine subjects
for iSubject = 1:nSubjects
    subjectID = subjectIDs{iSubject};
    fprintf('%s\n', subjectID)
    
    % load data
    subject = sprintf('%s_a1_tc100_soa1000-1250', subjectID);
    dataDir = pathToExpt('data');
    dataDir = sprintf('%s/%s/%s', dataDir, expName, subject(1:2));
    dataFile = dir(sprintf('%s/%s_run%02d*', dataDir, subject, run));
    load(sprintf('%s/%s', dataDir, dataFile(1).name))
    
    % add subject data to super subject structure
    SS.expt.trials = [SS.expt.trials; expt.trials];
    SS.expt.targetRotations = [SS.expt.targetRotations; expt.targetRotations];
    
    for iV = 1:3
        for iT = 1:2
            SS.results.totals.all{iV,iT} = [SS.results.totals.all{iV,iT}; results.totals.all{iV,iT}];
        end
    end
    
end

% can set this just once
SS.expt.p = expt.p;
SS.expt.trials_headers = expt.trials_headers;

%% save expt and results as for a normal subject
clear expt results
subject = 'SS_a1_tc100_soa1000-1250';
dataDir = pathToExpt('data');
dataDir = sprintf('%s/%s/%s', dataDir, expName, subject(1:2));
fileName = sprintf('%s/%s_run%02d_%s', dataDir, subject, run, datestr(now,'yyyymmdd'));

expt = SS.expt;
results = SS.results;
SS.subjectIDs = subjectIDs;

mkdir(dataDir)
save(fileName, 'expt', 'results', 'subject', 'run', 'SS')
