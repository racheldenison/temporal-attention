% rd_resampleTemporalAttention3TargetsGroup.m

% modified from rd_analyzeTemporalAttentionGroup.m

exptName = 'cbD15';
subjectInits = {'gb','xw','yz','jg','rd','ht','gb2','ds','ik','jp'};

run = 1;

nSamples = 1000;
% shuffleLabels = {'targetContrast','respInterval','cueValidity'};
shuffleLabels = {'cueValidity'};

nSubjects = numel(subjectInits);

% ibInfo = load('data/E5_nInvalidTrials_N10_20160822.mat');

%% Resample
for iSample = 1:nSamples
    fprintf('%d\n', iSample)
    
    %% Get data
    for iSubject = 1:nSubjects
        subjectInit = subjectInits{iSubject};
        
        dataDir = sprintf('%s/E5_T3_cbD15/%s', pathToExpt('data'), subjectInit);
        subjectID = sprintf('%s_%s', subjectInit, exptName);

        % load data from a given soa
        dataFile = dir(sprintf('%s/%s*_run%02d*_TemporalAttention3Targets*', ...
            dataDir, subjectID, run));
        if numel(dataFile)~=1
            sprintf('%s/%s*_run%02d*_TemporalAttention3Targets*', dataDir, subjectID, run)
            error('more or fewer than one matching data file')
        else
            load(sprintf('%s/%s', dataDir, dataFile.name))
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% shuffle labels and reanalyze %%%
        idx = zeros(size(shuffleLabels));
        for iLabel = 1:numel(shuffleLabels)
            idx(iLabel) = find(strcmp(expt.trials_headers, shuffleLabels{iLabel}));
        end
        newOrder = randperm(size(expt.trials,1));
        expt.trials(:,idx) = expt.trials(newOrder,idx);
        [expt results] = rd_analyzeTemporalAttention3Targets(expt);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        nT = numel(expt.p.respInterval);
        
        % read out the accuracy and rt
        % average across contrasts
        for iT = 1:nT 
            accDataC{iT}(:,iSubject) = mean(results.accMean{iT},2);
            rtDataC{iT}(:,iSubject) = mean(results.rtMean{iT},2);
        end
        
        % average across T1/T2/T3
        accDataCT = (accDataC{1} + accDataC{2} + accDataC{3})./2;
        rtDataCT = (rtDataC{1} + rtDataC{2} + rtDataC{3})./2;
        
        % break down invalid condition by which target was cued
        % (randomly generate cued target assignments for resampled data)
        invalid = rd_breakdownInvalidT3Resample(expt,results);
        nInvalidTrials(:,:,iSubject) = invalid.nTrials;
        ibNames = {'valid','neutral','invalid1','invalid2'};
        for iT = 1:nT
            ad = accDataC{iT}(:,iSubject);
            rd = rtDataC{iT}(:,iSubject);
            accDataIB{iT}(:,iSubject) = [ad([1 3]); invalid.accMean{iT,1}; invalid.accMean{iT,2}];
            rtDataIB{iT}(:,iSubject) = [rd([1 3]); invalid.rtMean{iT,1}; invalid.rtMean{iT,2}];
        end
    end
    
    % pairwise differences across subjects
    for iT = 1:nT 
        accDataCPairwise(1,iT,iSample) = mean(accDataC{iT}(1,:) - accDataC{iT}(2,:)); % VI
        accDataCPairwise(2,iT,iSample) = mean(accDataC{iT}(1,:) - accDataC{iT}(3,:)); % VN
        accDataCPairwise(3,iT,iSample) = mean(accDataC{iT}(3,:) - accDataC{iT}(2,:)); % NI
    end
    
    accDataCTPairwise(1,iSample) = mean(accDataCT(1,:) - accDataCT(2,:));
    accDataCTPairwise(2,iSample) = mean(accDataCT(1,:) - accDataCT(3,:));
    accDataCTPairwise(3,iSample) = mean(accDataCT(3,:) - accDataCT(2,:));
    
    for iT = 1:nT 
        rtDataCPairwise(1,iT,iSample) = mean(rtDataC{iT}(1,:) - rtDataC{iT}(2,:));
        rtDataCPairwise(2,iT,iSample) = mean(rtDataC{iT}(1,:) - rtDataC{iT}(3,:));
        rtDataCPairwise(3,iT,iSample) = mean(rtDataC{iT}(3,:) - rtDataC{iT}(2,:));
    end
    
    rtDataCTPairwise(1,iSample) = mean(rtDataCT(1,:) - rtDataCT(2,:));
    rtDataCTPairwise(2,iSample) = mean(rtDataCT(1,:) - rtDataCT(3,:));
    rtDataCTPairwise(3,iSample) = mean(rtDataCT(3,:) - rtDataCT(2,:));
    
    % pairwise differences  - invalid breakdown
    for iT = 1:nT 
        accDataIBPairwise(1,iT,iSample) = mean(accDataIB{iT}(1,:) - accDataIB{iT}(3,:)); % VI1
        accDataIBPairwise(2,iT,iSample) = mean(accDataIB{iT}(1,:) - accDataIB{iT}(4,:)); % VI2
        accDataIBPairwise(3,iT,iSample) = mean(accDataIB{iT}(1,:) - accDataIB{iT}(2,:)); % VN
        accDataIBPairwise(4,iT,iSample) = mean(accDataIB{iT}(2,:) - accDataIB{iT}(3,:)); % NI1
        accDataIBPairwise(5,iT,iSample) = mean(accDataIB{iT}(2,:) - accDataIB{iT}(4,:)); % NI2
    end

    for iT = 1:nT 
        rtDataIBPairwise(1,iT,iSample) = mean(rtDataIB{iT}(1,:) - rtDataIB{iT}(3,:)); % VI1
        rtDataIBPairwise(2,iT,iSample) = mean(rtDataIB{iT}(1,:) - rtDataIB{iT}(4,:)); % VI2
        rtDataIBPairwise(3,iT,iSample) = mean(rtDataIB{iT}(1,:) - rtDataIB{iT}(2,:)); % VN
        rtDataIBPairwise(4,iT,iSample) = mean(rtDataIB{iT}(2,:) - rtDataIB{iT}(3,:)); % NI1
        rtDataIBPairwise(5,iT,iSample) = mean(rtDataIB{iT}(2,:) - rtDataIB{iT}(4,:)); % NI2
    end
end

%% figures
vcLabels = {'VI','VN','NI'};
% acc
figure
for iVC = 1:3
    subplot(1,3,iVC)
    hist(squeeze(accDataCPairwise(iVC,:,:))')
    title(vcLabels{iVC})
end
rd_supertitle('acc');
rd_raiseAxis(gca);

figure
for iVC = 1:3
    subplot(1,3,iVC)
    hist(accDataCTPairwise(iVC,:))
    title(vcLabels{iVC})
end
rd_supertitle('acc');
rd_raiseAxis(gca);

figure
for iVC = 1:3
    subplot(1,3,iVC)
    hist(squeeze(rtDataCPairwise(iVC,:,:))')
    title(vcLabels{iVC})
end
rd_supertitle('RT');
rd_raiseAxis(gca);

figure
for iVC = 1:3
    subplot(1,3,iVC)
    hist(rtDataCTPairwise(iVC,:))
    title(vcLabels{iVC})
end
rd_supertitle('RT');
rd_raiseAxis(gca);

vcLabelsIB = {'VI1','VI2','VN','NI1','NI2'};
% acc
figure
for iVC = 1:5
    subplot(1,5,iVC)
    hist(squeeze(accDataIBPairwise(iVC,:,:))')
    title(vcLabelsIB{iVC})
end
rd_supertitle('acc');
rd_raiseAxis(gca);

% rt
figure
for iVC = 1:5
    subplot(1,5,iVC)
    hist(squeeze(rtDataIBPairwise(iVC,:,:))')
    title(vcLabelsIB{iVC})
end
rd_supertitle('RT');
rd_raiseAxis(gca);
