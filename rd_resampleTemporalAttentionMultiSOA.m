% rd_resampleTemporalAttentionMultiSOA.m

% modified from rd_resampleTemporalAttentionGroup.m

exptName = 'cbD6';
subjectInits = {'rd','hl','ho','vp'};
contrast = '*'; 
soa1 = [1000 1000 1000 1000 1000 1000 1000 1000 1000 1000];
soa2 = [1100 1150 1200 1250 1300 1350 1400 1450 1500 1800];
nSOAs = numel(soa1);

run = 98;

nSamples = 1000;
shuffleLabels = {'cueValidity'};

nSubjects = numel(subjectInits);

%% Resample
for iSubject = 1:nSubjects
    subjectInit = subjectInits{iSubject};
    fprintf('%s\n', subjectInit)
    
    for iSOA = 1:nSOAs
        fprintf('SOA %d\n', soa2(iSOA)-soa1(iSOA))
        exptStr = sprintf('%s*_*%s_soa%d-%d*', ...
            exptName, contrast, soa1(iSOA), soa2(iSOA));
        
        dataDir = sprintf('%s/E2_SOA_cbD6/%s', pathToExpt('data'), subjectInit(1:2));
        subjectID = sprintf('%s*_%s', subjectInit, exptStr);
        
        % load data from a given soa
        dataFile = dir(sprintf('%s/%s_run%02d*', ...
            dataDir, subjectID, run));
        if numel(dataFile)~=1
            sprintf('%s/%s_run%02d*', dataDir, subjectID, run)
            error('more or fewer than one matching data file')
        else
            load(sprintf('%s/%s', dataDir, dataFile.name))
        end
        
        for iSample = 1:nSamples
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% shuffle labels and reanalyze %%%
            idx = zeros(size(shuffleLabels));
            for iLabel = 1:numel(shuffleLabels)
                idx(iLabel) = find(strcmp(expt.trials_headers, shuffleLabels{iLabel}));
            end
            newOrder = randperm(size(expt.trials,1));
            expt.trials(:,idx) = expt.trials(newOrder,idx);
            [expt results] = rd_analyzeTemporalAttention(expt);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % read out the accuracy and rt
            % average across contrasts
            for iEL = 1:2 % early/late
                accData{iEL}(:,iSample,iSOA,iSubject) = results.accMean{iEL};
                rtData{iEL}(:,iSample,iSOA,iSubject) = results.rtMean{iEL};
            end
        end
    end
end

%% pairwise differences
for iEL = 1:2
    accDataPairwise(1,iEL,:,:,:) = accData{iEL}(1,:,:,:) - accData{iEL}(2,:,:,:); % VI
    accDataPairwise(2,iEL,:,:,:) = accData{iEL}(1,:,:,:) - accData{iEL}(3,:,:,:); % VN
    accDataPairwise(3,iEL,:,:,:) = accData{iEL}(3,:,:,:) - accData{iEL}(2,:,:,:); % NI
    
    rtDataPairwise(1,iEL,:,:,:) = rtData{iEL}(1,:,:,:) - rtData{iEL}(2,:,:,:);
    rtDataPairwise(2,iEL,:,:,:) = rtData{iEL}(1,:,:,:) - rtData{iEL}(3,:,:,:);
    rtDataPairwise(3,iEL,:,:,:) = rtData{iEL}(3,:,:,:) - rtData{iEL}(2,:,:,:);
end

accGroupPairwise = mean(accDataPairwise,5);
rtGroupPairwise = mean(rtDataPairwise,5);

%% save
% save data/E2_SOA_cbD6_randomizationTest_workspace_run98_N4_20160125.mat
