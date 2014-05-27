% rd_analyzeTemporalAttentionGroup.m

exptName = 'cb';
subjectInits = {'rd','ld','id'};
tilt = '*';
contrast = '64'; % '64-1*'
contrastIdx = 1; % only plot one contrast at a time
soa1 = 1000;
soa2 = 2500;

run = 9;

saveFigs = 1;

nSubjects = numel(subjectInits);
exptStr = sprintf('%s_tilt%s_tc%s_soa%d-%d*', ...
    exptName, tilt, contrast, soa1, soa2);

dataDir = pathToExpt('data');
% dataDir = [pathToExpt('data') '/pilot/rd'];

%% Get data
for iSubject = 1:nSubjects 
    subjectInit = subjectInits{iSubject};
    
%     dataDir = sprintf('%s/pilot/%s', pathToExpt('data'), subjectInit(1:2));
    subjectID = sprintf('%s_%s', subjectInit, exptStr);
    
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
        accData{iEL}(:,iSubject) = results.accMean{iEL}(:,contrastIdx);
        rtData{iEL}(:,iSubject) = results.rtMean{iEL}(:,contrastIdx);
    end
    
    % summary across subjects
    for iEL = 1:2 % early/late
        accMean(:,iEL) = mean(accData{iEL},2);
        accSte(:,iEL) = std(accData{iEL},0,2)./sqrt(nSubjects);
        rtMean(:,iEL) = mean(rtData{iEL},2);
        rtSte(:,iEL) = std(rtData{iEL},0,2)./sqrt(nSubjects);
    end 
end

p = expt.p;

%% Plot figs
intervalNames = {'early','late'};
accLims = [0.2 1];
rtLims = [0.3 1.3]; % [0.3 1.6]
xlims = [0 nSubjects+1];
colors = get(0,'DefaultAxesColorOrder');
axTitle = '';
[y, idx] = sort(p.cueValidity,2,'descend');

fig(1) = figure;
for iRI = 1:numel(p.respInterval)
    subplot(1,numel(p.respInterval),iRI);
    hold on
    plot(xlims, [0.5 0.5], '--k');
    
    p1 = bar(repmat((1:nSubjects)',1,numel(p.cueValidity)),...
        accData{iRI}(idx,:)');
%     colormap(colors(idx,:))
    colormap(flag(numel(p.cueValidity)));

    set(gca,'XTick',1:nSubjects)
    xlabel('subject')
    ylabel('acc')
    legend(p1, num2str(p.cueValidity(idx)'),'location','best')
    title(intervalNames{iRI})
    xlim(xlims)
    ylim(accLims)
    rd_supertitle(sprintf('%s run %d, N=%d', exptStr, run, nSubjects));
    rd_raiseAxis(gca);
end

fig(2) = figure;
for iRI = 1:numel(p.respInterval)
    subplot(1,numel(p.respInterval),iRI);
    
    p1 = bar(repmat((1:nSubjects)',1,numel(p.cueValidity)),...
        rtData{iRI}(idx,:)');
%     colormap(colors(idx,:))
    colormap(flag(numel(p.cueValidity)));
    
    set(gca,'XTick',1:nSubjects)
    xlabel('subject')
    ylabel('rt')
    legend(num2str(p.cueValidity(idx)'),'location','best')
    title(intervalNames{iRI})
    xlim(xlims)
    ylim(rtLims)
    box off
    rd_supertitle(sprintf('%s run %d, N=%d', exptStr, run, nSubjects));
    rd_raiseAxis(gca);
end

fig(3) = figure;
for iRI = 1:numel(p.respInterval)
    subplot(1,numel(p.respInterval),iRI);
    hold on
    
    plot([0 numel(p.cueValidity)+1], [0.5 0.5], '--k');
    b1 = bar(1:nSubjects, accMean(idx,iRI),'FaceColor',[.5 .5 .5]);
    p1 = errorbar(1:nSubjects, accMean(idx,iRI)', accSte(idx,iRI)','k','LineStyle','none');
    
    set(gca,'XTick',1:numel(p.cueValidity))
    set(gca,'XTickLabel', p.cueValidity(idx))
    xlabel('cue validity')
    ylabel('acc')
    title(intervalNames{iRI})
    ylim([0.3 max(accMean(:))*1.1])
    box off
    rd_supertitle(sprintf('%s run %d, N=%d', exptStr, run, nSubjects));
    rd_raiseAxis(gca);
end

fig(4) = figure;
for iRI = 1:numel(p.respInterval)
    subplot(1,numel(p.respInterval),iRI);
    hold on
    
    b1 = bar(1:nSubjects, rtMean(idx,iRI),'FaceColor',[.5 .5 .5]);
    p1 = errorbar(1:nSubjects, rtMean(idx,iRI)', rtSte(idx,iRI)','k','LineStyle','none');
    
    set(gca,'XTick',1:numel(p.cueValidity))
    set(gca,'XTickLabel', p.cueValidity(idx))
    xlabel('cue validity')
    ylabel('rt')
    title(intervalNames{iRI})
    ylim([min(rtMean(:))*0.9 max(rtMean(:))*1.1])
    box off
    rd_supertitle(sprintf('%s run %d, N=%d', exptStr, run, nSubjects));
    rd_raiseAxis(gca);
end

if saveFigs
    figNames = {'indivAcc','indivRT','groupAcc','groupRT'};
    rd_saveAllFigs(fig, figNames, sprintf('groupData_%s_run%02d_TemporalAttention_T1T2all_N%d', exptStr, run, nSubjects))
end
