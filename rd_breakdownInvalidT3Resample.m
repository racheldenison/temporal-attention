function invalid = rd_breakdownInvalidT3Resample(expt, results)

verbose = 0;

cuedIdx = strcmp(expt.trials_headers,'cuedInterval');
invalidIdx = expt.p.cueValidity==-1;
rtIdx = strcmp(expt.trials_headers,'rt');
correctIdx = strcmp(expt.trials_headers,'correct');

totalsAll = results.totals.all; % {iCV,iT}
nT = size(totalsAll,2);

% randomly generate cued target assignments (for resampled data)
for iT = 1:nT
    cued0 = totalsAll{invalidIdx,iT}(:,cuedIdx);
    cuedTargetOptions = setdiff(1:nT,iT);
    cued1 = repmat(cuedTargetOptions,1,round(numel(cued0)/numel(cuedTargetOptions)));
    cued{iT} = cued1(randperm(numel(cued0)));
end

if verbose
    fprintf('\n%s\n', expt.subjectID)
end
for iT = 1:nT
    if verbose
        fprintf('\nT%d invalid\n', iT)
    end
    c = cued{iT};
    cues = unique(c);
    for iC = 1:numel(cues)
        cIdx = c==cues(iC);
        if verbose
            fprintf('cue T%d: %d\n', cues(iC), nnz(cIdx))
        end
        
        invalid.cued(iT,iC) = cues(iC);
        invalid.nTrials(iT,iC) = nnz(cIdx);
        invalid.accMean{iT,iC} = nanmean(totalsAll{invalidIdx,iT}(cIdx,correctIdx));
        invalid.rtMean{iT,iC} = nanmean(totalsAll{invalidIdx,iT}(cIdx,rtIdx));
    end
end
    

