% rd_analyzeTemporalAttentionGroup.m

exptName = 'cb';
subjectInits = {'rd','ld','id','bl','ad','vp','ma','ty'};
tilt = '*';
contrast = '*'; % '64'; % 
contrastIdx = 1; % only plot one contrast at a time
soa1 = 1000;
soa2 = 1250;

run = 9;

saveFigs = 0;

nSubjects = numel(subjectInits);
exptStr = sprintf('%s_*%s_soa%d-%d*', ...
    exptName, contrast, soa1, soa2);

% dataDir = pathToExpt('data');
% dataDir = [pathToExpt('data') '/pilot/rd'];

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
    
    % read out the accuracy and rt
    for iEL = 1:2 % early/late
        accData{iEL}(:,iSubject) = results.accMean{iEL}(:,contrastIdx);
        rtData{iEL}(:,iSubject) = results.rtMean{iEL}(:,contrastIdx);
    end
    
    % also gather all the data for all contrasts
    for iContrast = 1:numel(expt.p.targetContrasts)
        for iEL = 1:2 % early/late
            accDataAll{iEL}(:,iSubject,iContrast) = results.accMean{iEL}(:,iContrast);
            rtDataAll{iEL}(:,iSubject,iContrast) = results.rtMean{iEL}(:,iContrast);
            
            % average across contrasts
            accDataC{iEL}(:,iSubject) = mean(results.accMean{iEL},2);
            rtDataC{iEL}(:,iSubject) = mean(results.rtMean{iEL},2);
        end
    end
end

% summary across subjects
for iEL = 1:2 % early/late
    accMean(:,iEL) = mean(accData{iEL},2);
    accSte(:,iEL) = std(accData{iEL},0,2)./sqrt(nSubjects);
    rtMean(:,iEL) = mean(rtData{iEL},2);
    rtSte(:,iEL) = std(rtData{iEL},0,2)./sqrt(nSubjects);
end

for iEL = 1:2 % early/late
    accMeanAll{iEL} = squeeze(mean(accDataAll{iEL},2));
    accSteAll{iEL} = squeeze(std(accDataAll{iEL},0,2)./sqrt(nSubjects));
    rtMeanAll{iEL} = squeeze(mean(rtDataAll{iEL},2));
    rtSteAll{iEL} = squeeze(std(rtDataAll{iEL},0,2)./sqrt(nSubjects));
end

p = expt.p;
nCV = numel(p.cueValidity);
tc = p.targetContrasts(contrastIdx)*100;

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
    
    p1 = bar(repmat((1:nSubjects)',1,nCV),...
        accData{iRI}(idx,:)');
%     colormap(colors(idx,:))
    colormap(flag(nCV));

    set(gca,'XTick',1:nSubjects)
    xlabel('subject')
    ylabel('acc')
    legend(p1, num2str(p.cueValidity(idx)'),'location','best')
    title(intervalNames{iRI})
    xlim(xlims)
    ylim(accLims)
    rd_supertitle(sprintf('%s run %d, contrast = %d%%, N=%d', exptStr, run, tc, nSubjects));
    rd_raiseAxis(gca);
end

fig(2) = figure;
for iRI = 1:numel(p.respInterval)
    subplot(1,numel(p.respInterval),iRI);
    
    p1 = bar(repmat((1:nSubjects)',1,nCV),...
        rtData{iRI}(idx,:)');
%     colormap(colors(idx,:))
    colormap(flag(nCV));
    
    set(gca,'XTick',1:nSubjects)
    xlabel('subject')
    ylabel('rt')
    legend(num2str(p.cueValidity(idx)'),'location','best')
    title(intervalNames{iRI})
    xlim(xlims)
    ylim(rtLims)
    box off
    rd_supertitle(sprintf('%s run %d, contrast = %d%%, N=%d', exptStr, run, tc, nSubjects));
    rd_raiseAxis(gca);
end

fig(3) = figure;
for iRI = 1:numel(p.respInterval)
    subplot(1,numel(p.respInterval),iRI);
    hold on
    
    plot([0 nCV+1], [0.5 0.5], '--k');
    b1 = bar(1:nCV, accMean(idx,iRI),'FaceColor',[.5 .5 .5]);
    p1 = errorbar(1:nCV, accMean(idx,iRI)', accSte(idx,iRI)','k','LineStyle','none');
    
    set(gca,'XTick',1:nCV)
    set(gca,'XTickLabel', p.cueValidity(idx))
    xlabel('cue validity')
    ylabel('acc')
    title(intervalNames{iRI})
    ylim([0.3 max(accMean(:))*1.1])
    ylim([.3 .9])
    box off
    rd_supertitle(sprintf('%s run %d, contrast = %d%%, N=%d', exptStr, run, tc, nSubjects));
    rd_raiseAxis(gca);
end

fig(4) = figure;
for iRI = 1:numel(p.respInterval)
    subplot(1,numel(p.respInterval),iRI);
    hold on
    
    b1 = bar(1:nCV, rtMean(idx,iRI),'FaceColor',[.5 .5 .5]);
    p1 = errorbar(1:nCV, rtMean(idx,iRI)', rtSte(idx,iRI)','k','LineStyle','none');
    
    set(gca,'XTick',1:nCV)
    set(gca,'XTickLabel', p.cueValidity(idx))
    xlabel('cue validity')
    ylabel('rt')
    title(intervalNames{iRI})
    ylim([min(rtMean(:))*0.9 max(rtMean(:))*1.1])
    box off
    rd_supertitle(sprintf('%s run %d, contrast = %d%%, N=%d', exptStr, run, tc, nSubjects));
    rd_raiseAxis(gca);
end

contrastLims = [p.targetContrasts(1)-0.05 p.targetContrasts(end)+0.05];
figure
for iRI = 1:numel(p.respInterval)
    subplot(1,numel(p.respInterval),iRI)
    hold on
    plot(contrastLims, [0.5 0.5], '--k');
    
    if numel(p.targetContrasts)>1
        p1 = errorbar(repmat(p.targetContrasts',1,numel(p.cueValidity)),...
            accMeanAll{iRI}', accSteAll{iRI}', '.', 'MarkerSize', 20);
    else
        for i = 1:length(accMean{iRI})
            p1(i) = errorbar(p.targetContrasts,...
                accMeanAll{iRI}(i), accSteAll{iRI}(i), '.', 'MarkerSize', 20);
            set(p1(i),'color', colors(i,:))
        end
    end
    xlabel('contrast')
    ylabel('acc')
    legend(p1, num2str(p.cueValidity'),'location','best')
    title(intervalNames{iRI})
    xlim(contrastLims)
    ylim(accLims)
    rd_supertitle(sprintf('%s run %d, contrast = %d%%, N=%d', exptStr, run, tc, nSubjects));
    rd_raiseAxis(gca);
    rd_supertitle(axTitle);
end

figure
rtLims = [0 2];
for iRI = 1:numel(p.respInterval)
    subplot(1,numel(p.respInterval),iRI)
    hold on
    
    if numel(p.targetContrasts)>1
        p1 = errorbar(repmat(p.targetContrasts',1,numel(p.cueValidity)),...
            rtMeanAll{iRI}', rtSteAll{iRI}', '.', 'MarkerSize', 20);
    else
        for i = 1:length(accMean{iRI})
            p1(i) = errorbar(p.targetContrasts,...
                rtMeanAll{iRI}(i), rtSteAll{iRI}(i), '.', 'MarkerSize', 20);
            set(p1(i),'color', colors(i,:))
        end
    end
    xlabel('contrast')
    ylabel('RT')
    legend(p1, num2str(p.cueValidity'),'location','best')
    title(intervalNames{iRI})
    xlim(contrastLims)
    ylim(rtLims)
    rd_supertitle(sprintf('%s run %d, contrast = %d%%, N=%d', exptStr, run, tc, nSubjects));
    rd_raiseAxis(gca);
    rd_supertitle(axTitle);
end

if saveFigs
    figNames = {'indivAcc','indivRT','groupAcc','groupRT'};
    rd_saveAllFigs(fig, figNames, sprintf('groupData_%s_run%02d_TemporalAttention_T1T2all_contrast%d_N%d', exptStr, run, tc, nSubjects))
end
