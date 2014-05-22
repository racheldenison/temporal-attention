% rd_plotTemporalAttentionMultiRun.m

subjectInit = 'rd';
exptName = 'cb';
tilt = '2';
contrast = '64';
contrastIdx = 1; % only plot one contrast at a time
soa1 = 1000;
soa2 = 1250;

runs = 1:3;

saveFigs = 1;

dataDir = pathToExpt('data');
% dataDir = [pathToExpt('data') '/pilot/rd'];

subjectID = sprintf('%s_%s_tilt%s_tc%s_soa%d-%d', ...
    subjectInit, exptName, tilt, contrast, soa1, soa2);

%% Get data
for iRun = 1:numel(runs)
    run = runs(iRun);
    
    % load data from a given soa
    dataFile = dir(sprintf('%s/%s_run%02d*', ...
        dataDir, subjectID, run));
    if numel(dataFile)~=1
        sprintf('%s/%s_run%02d*', dataDir, subjectID, run)
        error('more or fewer than one matching data file')
    else
        load(sprintf('%s/%s', dataDir, dataFile.name))
    end
    
    % read out the accuracy and rt
    for iEL = 1:2 % early/late
        accMean{iEL}(:,iRun) = results.accMean{iEL}(:,contrastIdx);
        accSte{iEL}(:,iRun) = results.accSte{iEL}(:,contrastIdx);
        rtMean{iEL}(:,iRun) = results.rtMean{iEL}(:,contrastIdx);
        rtSte{iEL}(:,iRun) = results.rtSte{iEL}(:,contrastIdx);
    end
end

p = expt.p;

%% Plot figs
intervalNames = {'early','late'};
accLims = [0.2 1];
rtLims = [0.3 1.6]; % [0.3 1.6]
xlims = [runs(1)-1 runs(end)+1];
colors = get(0,'DefaultAxesColorOrder');
axTitle = '';

fig(1) = figure;
for iRI = 1:numel(p.respInterval)
    subplot(1,numel(p.respInterval),iRI)
    hold on
    plot(xlims, [0.5 0.5], '--k');
    
    p1 = errorbar(repmat(runs',1,numel(p.cueValidity)),...
        accMean{iRI}', accSte{iRI}', '.-', 'MarkerSize', 20);

    set(gca,'XTick',runs)
    xlabel('run')
    ylabel('acc')
    legend(p1, num2str(p.cueValidity'),'location','best')
    title(intervalNames{iRI})
    xlim(xlims)
    ylim(accLims)
    rd_supertitle(subjectID);
    rd_raiseAxis(gca);
    rd_supertitle(axTitle);
end

fig(2) = figure;
for iRI = 1:numel(p.respInterval)
    subplot(1,numel(p.respInterval),iRI)
    
    p1 = errorbar(repmat(runs',1,numel(p.cueValidity)),...
        rtMean{iRI}', rtSte{iRI}', '.-', 'MarkerSize', 20);
    
    set(gca,'XTick',runs)
    xlabel('run')
    ylabel('rt')
    legend(num2str(p.cueValidity'),'location','best')
    title(intervalNames{iRI})
    xlim(xlims)
    ylim(rtLims)
    box off
    rd_supertitle(subjectID);
    rd_raiseAxis(gca);
    rd_supertitle(axTitle);
end

if saveFigs
    figNames = {'acc','rt'};
    rd_saveAllFigs(fig, figNames, sprintf('%s_byRun_TemporalAttention_T1T2all', subjectID))
end
