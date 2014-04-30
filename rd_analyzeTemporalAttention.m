function results = rd_analyzeTemporalAttention(expt, saveData, saveFigs)

if nargin < 3
    saveFigs = 0;
end
if nargin < 2
    saveData = 0;
end

%% Read out variables from expt
subjectID = expt.subjectID;
p = expt.p;
timing = expt.timing;
trials_headers = expt.trials_headers;
trials = expt.trials;
targetRotations = expt.targetRotations;

%% Get column idxs from trials_headers
targetContrastIdx = strcmp(trials_headers,'targetContrast');
respIntervalIdx = strcmp(trials_headers,'respInterval');
cueValidityIdx = strcmp(trials_headers,'cueValidity');
rtIdx = strcmp(trials_headers,'rt');
correctIdx = strcmp(trials_headers,'correct');

%% Analyze data
for iRI = 1:numel(p.respInterval)
    for iCV = 1:numel(p.cueValidity)
        for iTC = 1:numel(p.targetContrasts)
%             w0 = abs(diff(targetRotations,1,2))<10;
%             w0 = abs(diff(targetRotations,1,2))<190 & abs(diff(targetRotations,1,2))>170;

%             w = w0 & trials(:,respIntervalIdx)==iRI & trials(:,cueValidityIdx)==iCV & trials(:,targetContrastIdx)==iTC;
            w = trials(:,respIntervalIdx)==iRI & trials(:,cueValidityIdx)==iCV & trials(:,targetContrastIdx)==iTC;
            
            totals.all{iCV,iRI}(:,:,iTC) = trials(w,:);
            
            totals.means{iRI}(iCV,:,iTC) = mean(totals.all{iCV,iRI}(:,:,iTC),1);
            totals.stds{iRI}(iCV,:,iTC) = std(totals.all{iCV,iRI}(:,:,iTC),0,1);
            totals.stes{iRI}(iCV,:,iTC) = totals.stds{iRI}(iCV,:,iTC)./sqrt(size(totals.all{iCV,iRI}(:,:,iTC),1));
        end
    end
end

for iRI = 1:numel(p.respInterval)
    accMean{iRI} = squeeze(totals.means{iRI}(:,correctIdx,:)); % [validity x contrast]
    accSte{iRI} = squeeze(totals.stes{iRI}(:,correctIdx,:));
    
    rtMean{iRI} = squeeze(totals.means{iRI}(:,rtIdx,:));
    rtSte{iRI} = squeeze(totals.stes{iRI}(:,rtIdx,:));
end

%% Store data
results.totals = totals;
results.accMean = accMean;
results.accSte = accSte;
results.rtMean = rtMean;
results.rtSte = rtSte;
results.whenSaved = datestr(now);

%% Save data
if saveData
    fileName = sprintf('data/%s_TemporalAttention_%s.mat', subjectID, datestr(now, 'yyyymmdd'));
    save(fileName, 'expt', 'results')
end

%% Plot figs
intervalNames = {'early','late'};
accLims = [0.2 1];
rtLims = [0.3 1.6];
contrastLims = [p.targetContrasts(1)-0.05 p.targetContrasts(end)+0.05];

fig(1) = figure;
for iRI = 1:numel(p.respInterval)
    subplot(1,numel(p.respInterval),iRI)
    hold on
    plot(contrastLims, [0.5 0.5], '--k');
    p1 = errorbar(repmat(p.targetContrasts',1,numel(p.cueValidity)),...
        accMean{iRI}', accSte{iRI}', '.', 'MarkerSize', 20);
    xlabel('contrast')
    ylabel('acc')
    legend(p1, num2str(p.cueValidity'),'location','best')
    title(intervalNames{iRI})
    xlim(contrastLims)
    ylim(accLims)
    rd_supertitle(subjectID);
    rd_raiseAxis(gca);
end

fig(2) = figure;
for iRI = 1:numel(p.respInterval)
    subplot(1,numel(p.respInterval),iRI)
    errorbar(repmat(p.targetContrasts',1,numel(p.cueValidity)),...
        rtMean{iRI}', rtSte{iRI}', '.', 'MarkerSize', 20)
    xlabel('contrast')
    ylabel('rt')
    legend(num2str(p.cueValidity'),'location','best')
    title(intervalNames{iRI})
    xlim(contrastLims)
    ylim(rtLims)
    box off
    rd_supertitle(subjectID);
    rd_raiseAxis(gca);
end

%% Save figs
if saveFigs
    figNames = {'acc','rt'};
    rd_saveAllFigs(fig, figNames, [subjectID '_TemporalAttention'])
end


