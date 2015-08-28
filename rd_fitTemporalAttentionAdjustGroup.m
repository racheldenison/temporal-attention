% rd_fitTemporalAttentionAdjustGroup.m

%% setup
% subjectIDs = {'bl','rd','id','ec','ld','en','sj','ml','ca','jl','ew','jx'};
subjectIDs = {'id','ec','ld','en','sj','ml','ca','jl','ew','jx'};
run = 9;
nSubjects = numel(subjectIDs);

saveIndivData = 1;

bootstraps = 0; % 25:100

analysis = 'fit'; % 'fit' or 'compare'

%% get data
for iSubject = 1:nSubjects
    subjectID = subjectIDs{iSubject};
    fprintf('\n\n[%s]\n%s', datestr(now), subjectID)
    
    if any(bootstraps)
        for iBoot = 1:numel(bootstraps)
            bootRun = bootstraps(iBoot);
            fprintf('\nbootstrap %d', bootRun)
            switch analysis
                case 'fit'
                    rd_fitTemporalAttentionAdjust(subjectID, run, saveIndivData, bootRun);
                case 'compare'
                    rd_compareFitsTemporalAttentionAdjust(subjectID, run, saveIndivData, bootRun);
            end
        end
    else
        switch analysis
            case 'fit'
                [indivResults(iSubject).fit indivResults(iSubject).errors] = ...
                    rd_fitTemporalAttentionAdjust(subjectID, run, saveIndivData);
            case 'compare'
                [indivResults(iSubject).fit indivResults(iSubject).errors] = ...
                    rd_compareFitsTemporalAttentionAdjust(subjectID, run, saveIndivData);   
        end
    end
end