% rd_resampleTemporalAttentionAdjustGroup.m

%% setup
subjectIDs = {'bl','rd','id','ec','ld','en','sj','ml','ca','jl','ew','jx'};
run = 9;
nSubjects = numel(subjectIDs);

nSamples = 1000;

%% run resampling
for iSample = 1:nSamples
    tic
    fprintf('\n%d', iSample)
    
    for iSubject = 1:nSubjects
        subjectID = subjectIDs{iSubject};
%         fprintf('\n\n[%s]\n%s', datestr(now), subjectID)

        fit = rd_resampleTemporalAttentionAdjust(subjectID, run);
        
        for iV = 1:3
            for iT = 1:2
                p(iV,iT,:,iSubject) = fit(iV,iT).mle;
            end
        end
    end
    
    % calculate pairwise differences, mean across subjects
    pd(1,:,:,iSample) = mean(p(2,:,:,:) - p(1,:,:,:),4); % VI
    pd(2,:,:,iSample) = mean(p(3,:,:,:) - p(1,:,:,:),4); % VN
    pd(3,:,:,iSample) = mean(p(2,:,:,:) - p(3,:,:,:),4); % NI
    toc
end

%% plots
vcNames = {'VI','VN','NI'};
for iP = 1:size(pd,3)
    for iT = 1:2
        figure
        for iVC = 1:3
            subplot(1,3,iVC)
            hist(pd(iVC,iT,iP,:))
            title(vcNames{iVC})
        end
        rd_supertitle(sprintf('p%d T%d',iP,iT));
        rd_raiseAxis(gca);
    end
end


