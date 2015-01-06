% rd_fitTemporalAttentionAdjustGroup.m

%% setup
subjectIDs = {'bl','rd','id','ec','ld'};
run = 9;
nSubjects = numel(subjectIDs);

saveFigs = 0;

%% get data
for iSubject = 1:nSubjects
    subjectID = subjectIDs{iSubject};
    fprintf('\n\n[%s]\n%s', datestr(now), subjectID)
    indivResults(iSubject).fit = ...
        rd_fitTemporalAttentionAdjust(subjectID, run);
end