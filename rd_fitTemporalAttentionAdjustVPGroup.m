% rd_fitTemporalAttentionAdjustVPGroup.m

%% setup
% subjectIDs = {'bl','rd','id','ec','ld','en','sj','ml','ca','jl','ew','jx'};
subjectIDs = {'bl'};
run = 9;
nSubjects = numel(subjectIDs);

saveIndivData = 1;

bootstraps = 101:150; % 1:100

%% get data
for iSubject = 1:nSubjects
    subjectID = subjectIDs{iSubject};
    fprintf('\n\n[%s]\n%s', datestr(now), subjectID)
    
    if any(bootstraps)
        for iBoot = 1:numel(bootstraps)
            bootRun = bootstraps(iBoot);
            fprintf('\nbootstrap %d', bootRun)
            rd_fitTemporalAttentionAdjustVP(subjectID, run, saveIndivData, bootRun);
        end
    else
    [indivResults(iSubject).fit indivResults(iSubject).errors] = ...
        rd_fitTemporalAttentionAdjustVP(subjectID, run, saveIndivData);
    end
end