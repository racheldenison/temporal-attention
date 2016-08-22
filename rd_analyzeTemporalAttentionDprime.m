function results = rd_analyzeTemporalAttentionDprime(expt, results)

%% setup
accIdx = strcmp(expt.trials_headers,'correct');
targetStateIdx = strcmp(expt.trials_headers,'respTargetState');
targetState = expt.trials(:,targetStateIdx);
targetStates = unique(targetState);

totalsAll = results.totals.all;

nT = numel(expt.p.respInterval);

%% results separated by target state (ccw/cw)
for iT = 1:nT
    for iCV = 1:3
        for iState = 1:numel(targetStates)
            state = targetStates(iState);
            w = totalsAll{iCV,iT}(:,targetStateIdx)==state;
            meanByState{iState,iT}(iCV,:) = mean(totalsAll{iCV,iT}(w,:));
        end
    end
end

%% acc data by state
for iT = 1:nT
    for iState = 1:numel(targetStates)
        accStateMean{iT}(:,iState) = meanByState{iState,iT}(:,accIdx);
    end
end

%% calculate dprime
for iT = 1:nT
    [dprime{iT}, criterion{iT}] = rd_dprime(accStateMean{iT}(:,1), 1-accStateMean{iT}(:,2), '2afc','adjust');
end

%% store dprime and criterion
results.totals.dprime = dprime;
results.totals.crit = criterion;

