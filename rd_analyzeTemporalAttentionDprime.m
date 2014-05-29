function results = rd_analyzeTemporalAttentionDprime(expt, results)

%% setup
accIdx = strcmp(expt.trials_headers,'correct');
targetStateIdx = strcmp(expt.trials_headers,'respTargetState');
targetState = expt.trials(:,targetStateIdx);
targetStates = unique(targetState);

totalsAll = results.totals.all;

%% results separated by target state (ccw/cw)
for iEL = 1:2
    for iCV = 1:3
        for iState = 1:numel(targetStates)
            state = targetStates(iState);
            w = totalsAll{iCV,iEL}(:,targetStateIdx)==state;
            meanByState{iState,iEL}(iCV,:) = mean(totalsAll{iCV,iEL}(w,:));
        end
    end
end

%% acc data by state
for iEL = 1:2
    for iState = 1:numel(targetStates)
        accStateMean{iEL}(:,iState) = meanByState{iState,iEL}(:,accIdx);
    end
end

%% calculate dprime
for iEL = 1:2
    [dprime{iEL}, criterion{iEL}] = rd_dprime(accStateMean{iEL}(:,1), 1-accStateMean{iEL}(:,2), '2afc','adjust');
end

%% store dprime and criterion
results.totals.dprime = dprime;
results.totals.crit = criterion;

