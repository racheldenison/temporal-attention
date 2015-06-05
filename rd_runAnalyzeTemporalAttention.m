% rd_runAnalyzeTemporalAttention.m

subject = 'sl_cbD10_tilt*_tc16-64_soa1000-1800';
% subject = 'maPilot_cb_tilt5_soa1000-1250';
run = 8;

subjectID = sprintf('%s_run%02d', subject, run);

saveData = 1;
saveFigs = 1;
plotTimingFigs = 0;
saveTimingFigs = 0;
cleanRT = 0;

expName = 'E4_contrast_cbD10'; % 'E0_cb'
% dataDir = 'data';
% figDir = 'figures';
dataDir = pathToExpt('data');
figDir = pathToExpt('figures');
dataDir = sprintf('%s/%s/%s', dataDir, expName, subject(1:2));
figDir = sprintf('%s/%s/%s', figDir, expName, subject(1:2));

% load data file
dataFile = dir(sprintf('%s/%s_T*', dataDir, subjectID));
if numel(dataFile)~=1
    fprintf('\n%s/%s*', dataDir, subjectID)
    error('more or fewer than 1 matching data file')
else
    load(sprintf('%s/%s', dataDir, dataFile.name))
end

for t1t2 = {'all'} %{'same','diff'} %  
    [expt results] = rd_analyzeTemporalAttention(expt, 0, 0, 0, 0, t1t2{1}, cleanRT);
    
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
        if cleanRT
            figNames = {'cleanRTHist','acc','rt'};
        else
            figNames = {'acc','rt'};
        end
        rd_saveAllFigs([], figNames, [subjectID '_TemporalAttention_T1T2' t1t2{1}], figDir)
        close all % close so as not to interfere with next set of figures
    end 
end