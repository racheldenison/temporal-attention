% rd_runAnalyzeTemporalAttention.m

subject = 'rd_cb_tilt2_tc64_soa1000-2500';
run = 9;

subjectID = sprintf('%s_run%02d', subject, run);

saveData = 0;
saveFigs = 0;
plotTimingFigs = 0;
saveTimingFigs = 0;
cleanRT = 1;

dataDir = pathToExpt('data');
% dataDir = 'data';
% dataDir = [pathToExpt('data') '/pilot/rd'];

figDir = pathToExpt('figures');

% load data file
dataFile = dir(sprintf('%s/%s_T*', dataDir, subjectID));
if numel(dataFile)~=1
    fprintf('/n%s/%s*', dataDir, subjectID)
    error('more or fewer than 1 matching data file')
else
    load(sprintf('%s/%s', dataDir, dataFile.name))
end

for t1t2 = {'all'} % {'same','diff'} %
    rd_analyzeTemporalAttention(expt, 0, 0, 0, 0, t1t2{1}, cleanRT);
    
    % change subjectID for saving if cleaning RT
    if cleanRT
        subjectID = [subjectID '_RTx'];
    end

    % save data
    % saving data and figs separately in order to save them into the mcq
    % directory, and not locally
    if saveData
        fileName = sprintf('%s/%s_TemporalAttention_T1T2%s_%s.mat', dataDir, subjectID, t1t2{1}, datestr(now, 'yyyymmdd'));
        save(fileName, 'expt', 'results')
    end
    
    % save figs
    if saveFigs
        figNames = {'acc','rt'};
        rd_saveAllFigs([], figNames, [subjectID '_TemporalAttention_T1T2' t1t2{1}], figDir)
        close all % close so as not to interfere with next set of figures
    end 
end