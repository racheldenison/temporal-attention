% rd_resampleTemporalAttentionTradeoffs.m

% expNames = {'e5'};
expNames = {'e0','e3','e5'};
nExp = numel(expNames);

nSamples = 1000;

for iSample = 1:nSamples
    fprintf('%d\n', iSample)
    [pairedBC, pairNames, stackedBC, stackedNames] = ...
        rd_plotTemporalAttentionTradeoffsMean2(iSample);
    
    for iExp = 1:nExp
        expName = expNames{iExp};
        if strcmp(expName,'e5')
            pairedBCSamples.(expName)(:,:,:,iSample) = pairedBC.(expName);
            pairedBCAveSamples.(expName)(:,:,iSample) = mean(pairedBC.(expName),3);    
        else
            pairedBCSamples.(expName)(:,:,iSample) = pairedBC.(expName);
        end
        stackedBCSamples.(expName)(:,:,:,iSample) = stackedBC.(expName);
    end
end

for iExp = 1:nExp
    expName = expNames{iExp};
    if strcmp(expName,'e5')
        pairedBCCI.(expName) = prctile(pairedBCSamples.(expName),[16 84],4); % [2.5 97.5], [16 84]
        pairedBCAveCI.(expName) = prctile(pairedBCAveSamples.(expName),[16 84],3);
    else
        pairedBCCI.(expName) = prctile(pairedBCSamples.(expName),[16 84],3);
    end
    stackedBCCI.(expName) = prctile(stackedBCSamples.(expName),[16 84],4);
end