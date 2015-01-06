% rd_fitTemporalAttentionAdjustGroup.m

%% setup
subjectIDs = {'bl','rd','id','ec','ld'};
run = 19;
nSubjects = numel(subjectIDs);

saveIndivData = 1;

%% get data
for iSubject = 1:nSubjects
    subjectID = subjectIDs{iSubject};
    fprintf('\n\n[%s]\n%s', datestr(now), subjectID)
    [indivResults(iSubject).fit indivResults(iSubject).errors] = ...
        rd_fitTemporalAttentionAdjust(subjectID, run, saveIndivData);
end