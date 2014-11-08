% rd_plotTemporalAttentionMultiSOAGroup.m

%% Setup
subjectInits = {'rd','vp','hl'};
nSubjects = numel(subjectInits);
groupStr = sprintf('N=%d', nSubjects);

normalizeData = 0;

%% Get indiv subject data
for iSubject = 1:nSubjects
    subjectInit = subjectInits{iSubject};
    [acc rt t1t2soa p dp crit eff] = rd_plotTemporalAttentionMultiSOA(subjectInit);
    for iEL = 1:numel(acc)
        accData{iEL}(:,:,iSubject) = acc{iEL};
        rtData{iEL}(:,:,iSubject) = rt{iEL};
        dpData{iEL}(:,:,iSubject) = dp{iEL};
        critData{iEL}(:,:,iSubject) = crit{iEL};
        effData{iEL}(:,:,iSubject) = eff{iEL};
    end
end

%% Normalize data
if normalizeData
    for iEL = 1:numel(accData)
        accData{iEL} = normalizeDC(accData{iEL});
        rtData{iEL} = normalizeDC(rtData{iEL});
        dpData{iEL} = normalizeDC(dpData{iEL});
        critData{iEL} = normalizeDC(critData{iEL});
        effData{iEL} = normalizeDC(effData{iEL});
    end
end

%% Cuing effect and average across cue validities
for iEL = 1:numel(accData)
    accDataCueEff{iEL} = squeeze(accData{iEL}(1,:,:) - accData{iEL}(2,:,:));
    accDataCueAve{iEL} = squeeze(mean(accData{iEL},1));
end

%% Group summary stats
for iEL = 1:numel(accData)
    accMean{iEL} = mean(accData{iEL},3);
    rtMean{iEL} = mean(rtData{iEL},3);
    dpMean{iEL} = mean(dpData{iEL},3);
    critMean{iEL} = mean(critData{iEL},3);
    effMean{iEL} = mean(effData{iEL},3);
    
    accSte{iEL} = std(accData{iEL},0,3)./sqrt(nSubjects);
    rtSte{iEL} = std(rtData{iEL},0,3)./sqrt(nSubjects);
    dpSte{iEL} = std(dpData{iEL},0,3)./sqrt(nSubjects);
    critSte{iEL} = std(critData{iEL},0,3)./sqrt(nSubjects);
    effSte{iEL} = std(effData{iEL},0,3)./sqrt(nSubjects);
end

%% Plot figs
intervalNames = {'early','late'};
accLims = [0.2 1];
rtLims = [0.2 1.6];
dpLims = [-0.5 2.7];
critLims = [-1 1];
accDiffLims = [-0.15 0.25];
soaLims = [t1t2soa(1)-100 t1t2soa(end)+100];
colors = get(0,'DefaultAxesColorOrder');
axTitle = '';

fig(1) = figure;
for iRI = 1:numel(p.respInterval)
    subplot(1,numel(p.respInterval),iRI)
    hold on
    plot(soaLims, [0.5 0.5], '--k');
    
    p1 = errorbar(repmat(t1t2soa',1,numel(p.cueValidity)),...
        accMean{iRI}', accSte{iRI}', '.-', 'MarkerSize', 20);

    xlabel('soa')
    ylabel('acc')
    legend(p1, num2str(p.cueValidity'),'location','best')
    title(intervalNames{iRI})
    xlim(soaLims)
    ylim(accLims)
    rd_supertitle(subjectInits);
    rd_raiseAxis(gca);
    rd_supertitle(axTitle);
end

fig(2) = figure;
for iRI = 1:numel(p.respInterval)
    subplot(1,numel(p.respInterval),iRI)
    
    p1 = errorbar(repmat(t1t2soa',1,numel(p.cueValidity)),...
        rtMean{iRI}', rtSte{iRI}', '.-', 'MarkerSize', 20);
    
    xlabel('soa')
    ylabel('rt')
    legend(num2str(p.cueValidity'),'location','best')
    title(intervalNames{iRI})
    xlim(soaLims)
    ylim(rtLims)
    box off
    rd_supertitle(subjectInits);
    rd_raiseAxis(gca);
    rd_supertitle(axTitle);
end

fig(3) = figure;
for iRI = 1:numel(p.respInterval)
    subplot(1,numel(p.respInterval),iRI)
    hold on
    plot(soaLims, [0 0], '--k');
    
    p1 = errorbar(repmat(t1t2soa',1,numel(p.cueValidity)),...
        dpMean{iRI}', dpSte{iRI}', '.-', 'MarkerSize', 20);

    xlabel('soa')
    ylabel('dprime')
    legend(p1, num2str(p.cueValidity'),'location','best')
    title(intervalNames{iRI})
    xlim(soaLims)
    ylim(dpLims)
    rd_supertitle(subjectInits);
    rd_raiseAxis(gca);
    rd_supertitle(axTitle);
end

fig(4) = figure;
for iRI = 1:numel(p.respInterval)
    subplot(1,numel(p.respInterval),iRI)
    hold on
    plot(soaLims, [0 0], '--k');
    
    p1 = errorbar(repmat(t1t2soa',1,numel(p.cueValidity)),...
        critMean{iRI}', critSte{iRI}', '.-', 'MarkerSize', 20);

    xlabel('soa')
    ylabel('crit')
    legend(p1, num2str(p.cueValidity'),'location','best')
    title(intervalNames{iRI})
    xlim(soaLims)
    ylim(critLims)
    rd_supertitle(subjectInits);
    rd_raiseAxis(gca);
    rd_supertitle(axTitle);
end

fig(5) = figure;
for iRI = 1:numel(p.respInterval)
    subplot(1,numel(p.respInterval),iRI)
    hold on
    
    p1 = errorbar(repmat(t1t2soa',1,numel(p.cueValidity)),...
        effMean{iRI}', effSte{iRI}', '.-', 'MarkerSize', 20);

    xlabel('soa')
    ylabel('efficiency (dprime/rt)')
    legend(p1, num2str(p.cueValidity'),'location','best')
    title(intervalNames{iRI})
    xlim(soaLims)
%     ylim(effLims)
    rd_supertitle(subjectInits);
    rd_raiseAxis(gca);
    rd_supertitle(axTitle);
end

fig(6) = figure;
for iRI = 1:numel(p.respInterval)
    subplot(1,numel(p.respInterval),iRI)
    hold on
    plot(soaLims, [0 0], '--k');
    
    plot(t1t2soa, accDataCueEff{iRI})
    plot(t1t2soa, mean(accDataCueEff{iRI},2),'k','LineWidth',2)
    
    xlabel('soa')
    ylabel('cuing effect (accuracy valid-invalid)')
    title(intervalNames{iRI})
    xlim(soaLims)
    ylim(accDiffLims)
    rd_supertitle(subjectInits);
    rd_raiseAxis(gca);
    rd_supertitle(axTitle);
end

fig(7) = figure;
for iRI = 1:numel(p.respInterval)
    subplot(1,numel(p.respInterval),iRI)
    hold on
    plot(soaLims, [0 0], '--k');
    
    plot(t1t2soa, accDataCueAve{iRI})
    plot(t1t2soa, mean(accDataCueAve{iRI},2),'k','LineWidth',2)
    
    xlabel('soa')
    ylabel('acc')
    title(intervalNames{iRI})
    xlim(soaLims)
    ylim(accLims)
    rd_supertitle(subjectInits);
    rd_raiseAxis(gca);
    rd_supertitle(axTitle);
end


