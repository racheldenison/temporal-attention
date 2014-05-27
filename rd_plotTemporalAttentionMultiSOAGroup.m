% rd_plotTemporalAttentionMultiSOAGroup.m

%% Setup
subjectInits = {'rd','ld','id'};
nSubjects = numel(subjectInits);
groupStr = sprintf('N=%d', nSubjects);

%% Get indiv subject data
for iSubject = 1:nSubjects
    subjectInit = subjectInits{iSubject};
    [acc rt t1t2soa p] = rd_plotTemporalAttentionMultiSOA(subjectInit);
    for iEL = 1:numel(acc)
        accData{iEL}(:,:,iSubject) = acc{iEL};
        rtData{iEL}(:,:,iSubject) = rt{iEL};
    end
end

%% Group summary stats
for iEL = 1:numel(accData)
    accMean{iEL} = mean(accData{iEL},3);
    rtMean{iEL} = mean(rtData{iEL},3);
    
    accSte{iEL} = std(accData{iEL},0,3)./sqrt(nSubjects);
    rtSte{iEL} = std(rtData{iEL},0,3)./sqrt(nSubjects);
end

%% Plot figs
intervalNames = {'early','late'};
accLims = [0.2 1];
rtLims = [0.2 1.6];
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
