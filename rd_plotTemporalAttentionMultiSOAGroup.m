% rd_plotTemporalAttentionMultiSOAGroup.m

%% Setup
% subjectInits = {'rd','vp','hl','ho'};
% subjectInits = {'dg','sl','mr','ly','pv','ek','gk','md','ax'}; % 'ek'
subjectInits = {'dg','ly','ek'};
nSubjects = numel(subjectInits);
groupStr = sprintf('N=%d', nSubjects);

expName = 'E4'; % 'E2','E4'

contrast = 64;
T1T2Axis = ''; % 'same','diff'
extraSelection = ''; % 'sameOrient','diffOrient'
cleanRT = 1;

analStr1 = sprintf('contrast%d', contrast);
switch extraSelection
    case 'sameOrient'
        analStr2 = 'SO';
    case 'diffOrient'
        analStr2 = 'DO';
    otherwise
        analStr2 = '';
end
switch T1T2Axis
    case 'same'
        analStr3 = 'SA';
    case 'diff'
        analStr3 = 'DA';
    otherwise
        analStr3 = '';
end
analStr = sprintf('%s_%s%s', analStr1, analStr2, analStr3);
if strcmp(analStr(end),'_')
    analStr = analStr(1:end-1);
end

if cleanRT
    rtStr = '_RTx';
else
    rtStr = '';
end

figPrefix = sprintf('g%s_N%d_trained800_%s%s', expName, nSubjects, analStr, rtStr);

normalizeData = 0;
saveFigs = 0;

%% Get indiv subject data
for iSubject = 1:nSubjects
    subjectInit = subjectInits{iSubject};
    [acc rt t1t2soa p dp eff] = ...
        rd_plotTemporalAttentionMultiSOA(subjectInit, contrast, T1T2Axis, extraSelection, cleanRT);
    for iEL = 1:numel(acc)
        accData{iEL}(:,:,iSubject) = acc{iEL};
        rtData{iEL}(:,:,iSubject) = rt{iEL};
        dpData{iEL}(:,:,iSubject) = dp{iEL};
        effData{iEL}(:,:,iSubject) = eff{iEL};
    end
end

%% Normalize data
if normalizeData
    for iEL = 1:numel(accData)
        accData{iEL} = normalizeDC(accData{iEL});
        rtData{iEL} = normalizeDC(rtData{iEL});
        dpData{iEL} = normalizeDC(dpData{iEL});
        effData{iEL} = normalizeDC(effData{iEL});
    end
end

%% Cuing effect and average across cue validities
for iEL = 1:numel(accData)
    accDataCueEff{iEL} = squeeze(accData{iEL}(1,:,:) - accData{iEL}(2,:,:));
    accDataCueAve{iEL} = squeeze(mean(accData{iEL},1));
end

%% Cuing effect valid/invalid vs. neutral
for iEL = 1:numel(accData)
    accDataCueEffN{iEL}(:,:,1) = squeeze(accData{iEL}(1,:,:) - accData{iEL}(3,:,:));
    accDataCueEffN{iEL}(:,:,2) = squeeze(accData{iEL}(3,:,:) - accData{iEL}(2,:,:));
end

%% Group summary stats
for iEL = 1:numel(accData)
    accMean{iEL} = mean(accData{iEL},3);
    rtMean{iEL} = mean(rtData{iEL},3);
    dpMean{iEL} = mean(dpData{iEL},3);
    effMean{iEL} = mean(effData{iEL},3);
    
    accSte{iEL} = std(accData{iEL},0,3)./sqrt(nSubjects);
    rtSte{iEL} = std(rtData{iEL},0,3)./sqrt(nSubjects);
    dpSte{iEL} = std(dpData{iEL},0,3)./sqrt(nSubjects);
    effSte{iEL} = std(effData{iEL},0,3)./sqrt(nSubjects);
end

%% Plot figs
intervalNames = {'T1','T2'};
cueNames = {'valid','invalid','neutral'};
accLims = [0.2 1];
rtLims = [0.2 1.6];
dpLims = [-0.5 2.7];
effLims = [0 3.5];
accDiffLims = [-0.15 0.25];
accCueEffLims = [-0.3 0.35];
soaLims = [t1t2soa(1)-100 t1t2soa(end)+100];
colors = get(0,'DefaultAxesColorOrder');
axTitle = '';

% subjectColHSV = rgb2hsv(varycolor(nSubjects));
% subjectColors = hsv2rgb([subjectColHSV(:,1:2) subjectColHSV(:,3)*0.8]);
subjectColors = varycolor(nSubjects)*.9;
% subjectColors = colors((1:nSubjects)+2,:);

fig(1) = figure;
for iRI = 1:numel(p.respInterval)
    subplot(1,numel(p.respInterval),iRI)
    hold on
    plot(soaLims, [0.5 0.5], '--k');

    p1 = plot(repmat(t1t2soa',1,numel(p.cueValidity)),...
        accMean{iRI}', '-', 'LineWidth', 1.5);
    e1 = errorbar(repmat(t1t2soa',1,numel(p.cueValidity)),...
        accMean{iRI}', accSte{iRI}', '.', 'MarkerSize', 20, ...
        'LineWidth', 1);

    xlabel('soa')
    ylabel('acc')
    title(intervalNames{iRI})
    xlim(soaLims)
    ylim(accLims)
    legend(p1, cueNames,'location','best')
    rd_supertitle(subjectInits);
    rd_raiseAxis(gca);
    rd_supertitle(axTitle);
end

fig(2) = figure;
for iRI = 1:numel(p.respInterval)
    subplot(1,numel(p.respInterval),iRI)
    hold on

    p1 = plot(repmat(t1t2soa',1,numel(p.cueValidity)),...
        rtMean{iRI}', '-', 'LineWidth', 1.5);
    e1 = errorbar(repmat(t1t2soa',1,numel(p.cueValidity)),...
        rtMean{iRI}', rtSte{iRI}', '.', 'MarkerSize', 20, ...
        'LineWidth', 1);
    
    xlabel('soa')
    ylabel('rt')
    title(intervalNames{iRI})
    xlim(soaLims)
    ylim(rtLims)
    legend(p1, cueNames,'location','best')
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
    
    p1 = plot(repmat(t1t2soa',1,numel(p.cueValidity)),...
        dpMean{iRI}', '-','LineWidth', 1.5);
    e1 = errorbar(repmat(t1t2soa',1,numel(p.cueValidity)),...
        dpMean{iRI}', dpSte{iRI}', '.', 'MarkerSize', 20, ...
        'LineWidth', 1);

    xlabel('soa')
    ylabel('dprime')
    title(intervalNames{iRI})
    xlim(soaLims)
    ylim(dpLims)
    legend(p1, cueNames,'location','best')
    rd_supertitle(subjectInits);
    rd_raiseAxis(gca);
    rd_supertitle(axTitle);
end

fig(4) = figure;
for iRI = 1:numel(p.respInterval)
    subplot(1,numel(p.respInterval),iRI)
    hold on
    
    p1 = plot(repmat(t1t2soa',1,numel(p.cueValidity)),...
        effMean{iRI}', '-', 'LineWidth', 1.5);
    e1 = errorbar(repmat(t1t2soa',1,numel(p.cueValidity)),...
        effMean{iRI}', effSte{iRI}', '.', 'MarkerSize', 20, ...
        'LineWidth', 1);
    set(p1(2),'LineWidth',1.5)

    xlabel('soa')
    ylabel('efficiency (dprime/rt)')
    title(intervalNames{iRI})
    xlim(soaLims)
    ylim(effLims)
    legend(p1, cueNames,'location','best')
    rd_supertitle(subjectInits);
    rd_raiseAxis(gca);
    rd_supertitle(axTitle);
end

fig(5) = figure;
for iRI = 1:numel(p.respInterval)
    subplot(1,numel(p.respInterval),iRI)
    set(gca,'ColorOrder',subjectColors)
    hold all
    plot(soaLims, [0 0], '--k');
    
    p1 = plot(t1t2soa, accDataCueEff{iRI},'LineWidth',1.5);
    plot(t1t2soa, mean(accDataCueEff{iRI},2),'k','LineWidth',2.5)
    
    xlabel('soa')
    ylabel('cuing effect (accuracy valid-invalid)')
    title(intervalNames{iRI})
    xlim(soaLims)
    ylim(accDiffLims)
    rd_supertitle(subjectInits);
    rd_raiseAxis(gca);
    rd_supertitle(axTitle);
end

fig(6) = figure;
for iRI = 1:numel(p.respInterval)
    subplot(1,numel(p.respInterval),iRI)
    set(gca,'ColorOrder',subjectColors)
    hold all
    plot(soaLims, [0 0], '--k');
    
    plot(t1t2soa, accDataCueAve{iRI},'LineWidth',1.5)
    plot(t1t2soa, mean(accDataCueAve{iRI},2),'k','LineWidth',2.5)
    
    xlabel('soa')
    ylabel('acc average')
    title(intervalNames{iRI})
    xlim(soaLims)
    ylim(accLims)
    rd_supertitle(subjectInits);
    rd_raiseAxis(gca);
    rd_supertitle(axTitle);
end


fig(7) = figure;
hold on
plot(soaLims, [0 0], '--k')
p1 = [];
p1(1) = plot(t1t2soa, mean(accDataCueEff{1},2),'-','LineWidth',2);
p1(2) = plot(t1t2soa, mean(accDataCueEff{2},2),'r-','LineWidth',2);

e1(1) = errorbar(t1t2soa, mean(accDataCueEff{1},2), ...
    std(accDataCueEff{1},0,2)./sqrt(nSubjects),...
    '.','MarkerSize', 20,'LineWidth', 1);
e1(2) = errorbar(t1t2soa, mean(accDataCueEff{2},2), ...
    std(accDataCueEff{1},0,2)./sqrt(nSubjects),...
    '.r','MarkerSize', 20,'LineWidth', 1);

plot(t1t2soa, mean(accDataCueEff{1},2),'o','LineWidth',1,'MarkerSize',8,'MarkerFaceColor','w')
plot(t1t2soa, mean(accDataCueEff{2},2),'rs','LineWidth',1,'MarkerSize',8,'MarkerFaceColor','w')

legend(p1, intervalNames)
xlabel('soa')
ylabel('cuing effect (accuracy valid-invalid)')
xlim(soaLims)
ylim(accCueEffLims)


fig(8) = figure;
hold on
plot(soaLims, [0.5 0.5], '--k');
p1 = [];
p1(1) = plot(t1t2soa, mean(accDataCueAve{1},2),'-','LineWidth',2);
p1(2) = plot(t1t2soa, mean(accDataCueAve{2},2),'r-','LineWidth',2);

e1(1) = errorbar(t1t2soa, mean(accDataCueAve{1},2), ...
    std(accDataCueAve{1},0,2)./sqrt(nSubjects),...
    '.','MarkerSize', 20,'LineWidth', 1);
e1(2) = errorbar(t1t2soa, mean(accDataCueAve{2},2), ...
    std(accDataCueAve{1},0,2)./sqrt(nSubjects),...
    '.r','MarkerSize', 20,'LineWidth', 1);

plot(t1t2soa, mean(accDataCueAve{1},2),'o','LineWidth',1,'MarkerSize',8,'MarkerFaceColor','w')
plot(t1t2soa, mean(accDataCueAve{2},2),'rs','LineWidth',1,'MarkerSize',8,'MarkerFaceColor','w')

legend(p1, intervalNames)
xlabel('soa')
ylabel('acc average')
xlim(soaLims)
ylim(accLims)

viNames = {'valid - neutral','neutral - invalid'};
fig(9) = figure;
for iVI = 1:2
    subplot(2,1,iVI)
    hold on
    plot(soaLims, [0 0], '--k')
    p1 = [];
    p1(1,:) = plot(t1t2soa, squeeze(mean(accDataCueEffN{1}(:,:,iVI),2)),'-','LineWidth',2);
    p1(2,:) = plot(t1t2soa, squeeze(mean(accDataCueEffN{2}(:,:,iVI),2)),'-r','LineWidth',2);
    
    e1 = [];
    e1(1,:) = errorbar(t1t2soa, squeeze(mean(accDataCueEffN{1}(:,:,iVI),2)), ...
        squeeze(std(accDataCueEffN{1}(:,:,iVI),0,2))./sqrt(nSubjects),...
        '.','MarkerSize', 20,'LineWidth', 1);
    e1(2,:) = errorbar(t1t2soa, squeeze(mean(accDataCueEffN{2}(:,:,iVI),2)), ...
        squeeze(std(accDataCueEffN{1}(:,:,iVI),0,2))./sqrt(nSubjects),...
        '.r','MarkerSize', 20,'LineWidth', 1);
    
    plot(t1t2soa, squeeze(mean(accDataCueEffN{1}(:,:,iVI),2)),'o','LineWidth',1,'MarkerSize',8,'MarkerFaceColor','w')
    plot(t1t2soa, squeeze(mean(accDataCueEffN{2}(:,:,iVI),2)),'rs','LineWidth',1,'MarkerSize',8,'MarkerFaceColor','w')
    
    legend(p1(:,1), intervalNames)
    xlabel('soa')
    ylabel('cuing effect (accuracy)')
    xlim(soaLims)
    title(viNames{iVI})
end

%% set figure properties
for iF = 1:numel(fig)
    % set font size of titles, axis labels, and legends
    % see also: http://stackoverflow.com/questions/8934468/changing-fonts-in-matlab-plots
    set(findall(fig(iF),'type','text'),'FontSize',14)
end

%% Save figs
% figPrefix = [figPrefix '_sameAxDiffOrient'];
if saveFigs
    figNames = {'acc','rt','dprime','eff','accCueEffect','accCueAve','accCueEffectOverlay','accCueAveOverlay','accCueEffectNOverlay'};
    rd_saveAllFigs(fig, figNames, figPrefix, [])
end
