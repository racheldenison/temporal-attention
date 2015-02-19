function nonTargetErrors = getAdjustNonTargetErrors(subjectID, run)
%
% function nonTargetErrors = getAdjustNonTargetErrors(subjectID, run)

[e, to, nto, tod, po, pod, r] = ...
    rd_plotTemporalAttentionAdjustErrors(subjectID, run, 0);

for iEL = 1:2
    for iV = 1:3
        responses = r{iV,iEL};
        nonTargetOrients = nto{iV,iEL};
        
        nonTargetErrors{iV,iEL} = ...
            calculateOrientationResponseError2(nonTargetOrients, responses);
    end
end