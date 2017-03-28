% rd_resampleTemporalAttentionGroupBC.m
tic
% modified from rd_analyzeTemporalAttentionGroup.m

exptName = 'cb';
subjectInits = {'rd','ld','id','bl','ad','vp','ma','ty','zw','ec'};
contrast = '*'; 
soa1 = 1000;
soa2 = 1250;

run = 9;

resampleOption = 'bootstrap'; % 'bootstrap','permutation'
nSamples = 1000;
% shuffleLabels = {'targetContrast','respInterval','cueValidity'};
shuffleLabels = {'cueValidity'};

nSubjects = numel(subjectInits);
exptStr = sprintf('%s_*%s_soa%d-%d*', ...
    exptName, contrast, soa1, soa2);

%% Resample
for iSample = 1:nSamples
    fprintf('%d\n', iSample)
    
    %% Get data
    for iSubject = 1:nSubjects
        subjectInit = subjectInits{iSubject};
        
        dataDir = sprintf('%s/E0_cb/%s', pathToExpt('data'), subjectInit(1:2));
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
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        switch resampleOption
            case 'bootstrap'
                %%% resample and reanalyze %%%
                cvIdx = find(strcmp(expt.trials_headers, 'cueValidity'));
                tIdx = find(strcmp(expt.trials_headers, 'respInterval'));
                for iCV = 1:3
                    wcv = expt.trials(:,cvIdx)==iCV;
                    for iT = 1:2
                        wt = expt.trials(:,tIdx)==iT;
                        w = wcv & wt;
                        vals = expt.trials(w,:);
                        resamp = RandSample(1:nnz(w),[1 nnz(w)]);
                        expt.trials(w,:) = vals(resamp,:);    
                    end
                end
                [expt results] = rd_analyzeTemporalAttention(expt);     
                
            case 'permutation'
                %%% shuffle labels and reanalyze %%%
                idx = zeros(size(shuffleLabels));
                for iLabel = 1:numel(shuffleLabels)
                    idx(iLabel) = find(strcmp(expt.trials_headers, shuffleLabels{iLabel}));
                end
                newOrder = randperm(size(expt.trials,1));
                expt.trials(:,idx) = expt.trials(newOrder,idx);
                [expt results] = rd_analyzeTemporalAttention(expt);
            otherwise
                error('resample option not recognized')
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % read out the accuracy and rt
        % average across contrasts
        for iEL = 1:2 % early/late
            accDataC{iEL}(:,iSubject) = mean(results.accMean{iEL},2);
            rtDataC{iEL}(:,iSubject) = mean(results.rtMean{iEL},2);
        end
        
        % average across T1/T2
        accDataCT = (accDataC{1} + accDataC{2})./2;
        rtDataCT = (rtDataC{1} + rtDataC{2})./2;
    end
    
    % group mean
    for iEL = 1:2 % early/late
        accMeanC(:,iEL,iSample) = mean(accDataC{iEL},2);
        rtMeanC(:,iEL,iSample) = mean(rtDataC{iEL},2);
    end
    
    % pairwise differences across subjects
    for iEL = 1:2 
        accDataCPairwise(1,iEL,iSample) = mean(accDataC{iEL}(1,:) - accDataC{iEL}(2,:)); % VI
        accDataCPairwise(2,iEL,iSample) = mean(accDataC{iEL}(1,:) - accDataC{iEL}(3,:)); % VN
        accDataCPairwise(3,iEL,iSample) = mean(accDataC{iEL}(3,:) - accDataC{iEL}(2,:)); % NI
    end
    
    accDataCTPairwise(1,iSample) = mean(accDataCT(1,:) - accDataCT(2,:));
    accDataCTPairwise(2,iSample) = mean(accDataCT(1,:) - accDataCT(3,:));
    accDataCTPairwise(3,iSample) = mean(accDataCT(3,:) - accDataCT(2,:));
    
    for iEL = 1:2
        rtDataCPairwise(1,iEL,iSample) = mean(rtDataC{iEL}(1,:) - rtDataC{iEL}(2,:));
        rtDataCPairwise(2,iEL,iSample) = mean(rtDataC{iEL}(1,:) - rtDataC{iEL}(3,:));
        rtDataCPairwise(3,iEL,iSample) = mean(rtDataC{iEL}(3,:) - rtDataC{iEL}(2,:));
    end
    
    rtDataCTPairwise(1,iSample) = mean(rtDataCT(1,:) - rtDataCT(2,:));
    rtDataCTPairwise(2,iSample) = mean(rtDataCT(1,:) - rtDataCT(3,:));
    rtDataCTPairwise(3,iSample) = mean(rtDataCT(3,:) - rtDataCT(2,:));
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

toc