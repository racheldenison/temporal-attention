% rd_analyzeTemporalAttentionWorkspace.m

%% load workspace

%% setup
% not completed trials
nc = trials(:,11)==0;
trials(nc,7:end) = NaN;

%% Get column idxs from trials_headers
targetContrastIdx = strcmp(trials_headers,'targetContrast');
respIntervalIdx = strcmp(trials_headers,'respInterval');
cueValidityIdx = strcmp(trials_headers,'cueValidity');
rtIdx = strcmp(trials_headers,'rt');
correctIdx = strcmp(trials_headers,'correct');

%% selection
wAx = ones(size(targetRotations,1),1);
axTitle = '';

%% standard error by trial
for iRI = 1:numel(p.respInterval)
    for iCV = 1:numel(p.cueValidity)
        for iTC = 1:numel(p.targetContrasts)
            if exist('w0','var')
                w = w0 & wAx & trials(:,respIntervalIdx)==iRI & trials(:,cueValidityIdx)==iCV & trials(:,targetContrastIdx)==iTC;
            else
                w = wAx & trials(:,respIntervalIdx)==iRI & trials(:,cueValidityIdx)==iCV & trials(:,targetContrastIdx)==iTC;
            end
            
            nt = nnz(~isnan(trials(w,end)));
            try
                totals.all{iCV,iRI}(:,:,iTC) = trials(w,:);
                
                totals.means{iRI}(iCV,:,iTC) = nanmean(totals.all{iCV,iRI}(:,:,iTC),1);
                totals.stds{iRI}(iCV,:,iTC) = nanstd(totals.all{iCV,iRI}(:,:,iTC),0,1);
                totals.stes{iRI}(iCV,:,iTC) = totals.stds{iRI}(iCV,:,iTC)./sqrt(nt);
            catch
                fprintf('\nwarning: unequal numbers of trials in different conditions\n\n')
                totals.all{iCV,iRI}(:,:,iTC) = NaN;
                
                totals.means{iRI}(iCV,:,iTC) = nanmean(trials(w,:),1);
                totals.stds{iRI}(iCV,:,iTC) = nanstd(trials(w,:),0,1);
                totals.stes{iRI}(iCV,:,iTC) = totals.stds{iRI}(iCV,:,iTC)./sqrt(nt);
            end
        end
    end
end

%% Acc and RT means
for iRI = 1:numel(p.respInterval)
    accMean{iRI} = squeeze(totals.means{iRI}(:,correctIdx,:)); % [validity x contrast]
    accSte{iRI} = squeeze(totals.stes{iRI}(:,correctIdx,:));
    
    rtMean{iRI} = squeeze(totals.means{iRI}(:,rtIdx,:));
    rtSte{iRI} = squeeze(totals.stes{iRI}(:,rtIdx,:));
end

%% Plot figs
intervalNames = {'T1','T2','T3'};
accLims = [0.2 1];
rtLims = [0 1]; %[0.2 1.6]; % [0.8 2.2];
contrastLims = [p.targetContrasts(1)-0.05 p.targetContrasts(end)+0.05];
colors = get(0,'DefaultAxesColorOrder');

figure
imagesc(trials)

fig(1) = figure;
for iRI = 1:numel(p.respInterval)
    subplot(1,numel(p.respInterval),iRI)
    hold on
    plot(contrastLims, [0.5 0.5], '--k');
    
    if numel(p.targetContrasts)>1
        p1 = errorbar(repmat(p.targetContrasts',1,numel(p.cueValidity)),...
            accMean{iRI}', accSte{iRI}', '.', 'MarkerSize', 20);
    else
        for i = 1:length(accMean{iRI})
            p1(i) = errorbar(p.targetContrasts,...
                accMean{iRI}(i), accSte{iRI}(i), '.', 'MarkerSize', 20);
            set(p1(i),'color', colors(i,:))
        end
    end
    xlabel('contrast')
    ylabel('acc')
    legend(p1, num2str(p.cueValidity'),'location','best')
    title(intervalNames{iRI})
    xlim(contrastLims)
    ylim(accLims)
    rd_supertitle(subjectID);
    rd_raiseAxis(gca);
    rd_supertitle(axTitle);
end

fig(2) = figure;
for iRI = 1:numel(p.respInterval)
    subplot(1,numel(p.respInterval),iRI)
    hold on
    
    if numel(p.targetContrasts)>1
        p1 = errorbar(repmat(p.targetContrasts',1,numel(p.cueValidity)),...
            rtMean{iRI}', rtSte{iRI}', '.', 'MarkerSize', 20);
    else
        for i = 1:length(rtMean{iRI})
            p1(i) = errorbar(p.targetContrasts,...
                rtMean{iRI}(i), rtSte{iRI}(i), '.', 'MarkerSize', 20);
            set(p1(i),'color', colors(i,:))
        end
    end

    xlabel('contrast')
    ylabel('rt')
    legend(num2str(p.cueValidity'),'location','best')
    title(intervalNames{iRI})
    xlim(contrastLims)
    ylim(rtLims)
    box off
    rd_supertitle(subjectID);
    rd_raiseAxis(gca);
    rd_supertitle(axTitle);
end
