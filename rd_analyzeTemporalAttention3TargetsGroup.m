% rd_analyzeTemporalAttention3TargetsGroup.m

exptName = 'cbD15';
subjectInits = {'gb','xw','yz','jg','rd','ht','gb2','ds','ik','jp'};
contrastIdx = 1;

run = 1;

doDprime = 0;
normalizeData = 0;
if normalizeData
    normStr = '_n';
else
    normStr = '';
end

saveFigs = 0;

nSubjects = numel(subjectInits);

% dataDir = pathToExpt('data');
% dataDir = [pathToExpt('data') '/pilot/rd'];

doRandomizationTests = 0; % requries having generated empirical null distributions

%% Get data
for iSubject = 1:nSubjects 
    subjectInit = subjectInits{iSubject};
    
    dataDir = sprintf('%s/E5_T3_cbD15/%s', pathToExpt('data'), subjectInit);
    subjectID = sprintf('%s_%s', subjectInit, exptName);
    
    % load data from a given soa
    dataFile = dir(sprintf('%s/%s*_run%02d*_TemporalAttention3Targets*', ...
        dataDir, subjectID, run));
    if numel(dataFile)~=1
        sprintf('%s/%s*_run%02d*_TemporalAttention3Targets*', dataDir, subjectID, run)
        error('more or fewer than one matching data file')
    else
        load(sprintf('%s/%s', dataDir, dataFile.name))
    end
    
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%% if you want to reanalyze, do it here %%%
%     T1T2Axis = 'same';
%     [expt results] = rd_analyzeTemporalAttention3Targets(expt, 0, 0, 0, 0, T1T2Axis, 0);
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    nT = numel(expt.p.respInterval);

    % read out the accuracy and rt
    for iT = 1:nT % early/late
        accData{iT}(:,iSubject) = results.accMean{iT}(:,contrastIdx);
        rtData{iT}(:,iSubject) = results.rtMean{iT}(:,contrastIdx);
        
        % replace accData with dprime if requested
        if doDprime
%             accData{iT}(:,iSubject) = ...
%                 rd_dprime(accData{iT}(:,iSubject),[],'2afc','adjust');
            [dprime, criterion] = rd_dprimeTemporalAttention(expt.trials, expt.trials_headers);
            accData{iT}(:,iSubject) = dprime(:,iT);
        end
    end
    
    % also gather all the data for all contrasts
    %     for iContrast = 1:numel(expt.p.targetContrasts)
    for iT = 1:nT % early/late
        %             accDataAll{iT}(:,iSubject,iContrast) = results.accMean{iT}(:,iContrast);
        %             rtDataAll{iT}(:,iSubject,iContrast) = results.rtMean{iT}(:,iContrast);
        
        % average across contrasts
        accDataC{iT}(:,iSubject) = mean(results.accMean{iT},2);
        rtDataC{iT}(:,iSubject) = mean(results.rtMean{iT},2);
        
        % replace accData with dprime if requested
        if doDprime
%             accDataC{iT}(:,iSubject) = ...
%                 rd_dprime(accDataC{iT}(:,iSubject),[],'2afc','adjust');
            [dprime, criterion] = rd_dprimeTemporalAttention(expt.trials, expt.trials_headers);
            accDataC{iT}(:,iSubject) = dprime(:,iT);
        end
    end
    %     end
    
    % break down invalid condition by which target was cued
    invalid = rd_breakdownInvalidT3(expt,results);
    nInvalidTrials(:,:,iSubject) = invalid.nTrials;
    ibNames = {'valid','neutral','invalid1','invalid2'};
    for iT = 1:nT
        ad = accData{iT}(:,iSubject);
        rd = rtData{iT}(:,iSubject);
        accDataIB{iT}(:,iSubject) = [ad([1 3]); invalid.accMean{iT,1}; invalid.accMean{iT,2}];
        rtDataIB{iT}(:,iSubject) = [rd([1 3]); invalid.rtMean{iT,1}; invalid.rtMean{iT,2}];
    end
end

% inverse efficiency
for iT = 1:nT
    invEffDataC{iT} = rtDataC{iT}./accDataC{iT};
end

%% Summary across subjects
for iT = 1:nT
    accMean(:,iT) = mean(accData{iT},2);
    accSte(:,iT) = std(accData{iT},0,2)./sqrt(nSubjects);
    rtMean(:,iT) = mean(rtData{iT},2);
    rtSte(:,iT) = std(rtData{iT},0,2)./sqrt(nSubjects);
end

% invalid breakdown
for iT = 1:nT 
    accMeanIB(:,iT) = mean(accDataIB{iT},2);
    accSteIB(:,iT) = std(accDataIB{iT},0,2)./sqrt(nSubjects);
    rtMeanIB(:,iT) = mean(rtDataIB{iT},2);
    rtSteIB(:,iT) = std(rtDataIB{iT},0,2)./sqrt(nSubjects);
end

%% Normalize
if normalizeData
    % use the Morey 2008 correction
    % avarage across contrast
    a = cat(3, accDataC{1}, accDataC{2}, accDataC{3});
    b = normalizeDC(shiftdim(a,2));
    c = shiftdim(b,1);
    accDataC{1} = c(:,:,1);
    accDataC{2} = c(:,:,2);
    accDataC{3} = c(:,:,3);
    
    a = cat(3, rtDataC{1}, rtDataC{2}, rtDataC{3});
    b = normalizeDC(shiftdim(a,2));
    c = shiftdim(b,1);
    rtDataC{1} = c(:,:,1);
    rtDataC{2} = c(:,:,2);
    rtDataC{3} = c(:,:,3);
    
    [fixed1 fixed2 N] = size(b);
    M = fixed1*fixed2;
    morey = M/(M-1);
    
    % invalid breakdown
    a = cat(3, accDataIB{1}, accDataIB{2}, accDataIB{3});
    b = normalizeDC(shiftdim(a,2));
    c = shiftdim(b,1);
    accDataIB{1} = c(:,:,1);
    accDataIB{2} = c(:,:,2);
    accDataIB{3} = c(:,:,3);
    
    a = cat(3, rtDataIB{1}, rtDataIB{2}, rtDataIB{3});
    b = normalizeDC(shiftdim(a,2));
    c = shiftdim(b,1);
    rtDataIB{1} = c(:,:,1);
    rtDataIB{2} = c(:,:,2);
    rtDataIB{3} = c(:,:,3);
    
    [fixed1 fixed2 N] = size(b);
    M = fixed1*fixed2;
    moreyIB = M/(M-1);
    
    for iT = 1:nT
        accMeanC(:,iT) = mean(accDataC{iT},2);
        accSteC(:,iT) = sqrt(morey*var(accDataC{iT},0,2)./(nSubjects));
        rtMeanC(:,iT) = mean(rtDataC{iT},2);
        rtSteC(:,iT) = sqrt(morey*var(rtDataC{iT},0,2)./(nSubjects));
        
        accMeanIB(:,iT) = mean(accDataIB{iT},2);
        accSteIB(:,iT) = sqrt(moreyIB*var(accDataIB{iT},0,2)./(nSubjects));
        rtMeanIB(:,iT) = mean(rtDataIB{iT},2);
        rtSteIB(:,iT) = sqrt(moreyIB*var(rtDataIB{iT},0,2)./(nSubjects));
    end
else
    for iT = 1:nT % early/late
        %     if normalizeData
        %         accDataC{iT} = normalizeDC(accDataC{iT});
        %         rtDataC{iT} = normalizeDC(rtDataC{iT});
        %     end
        accMeanC(:,iT) = mean(accDataC{iT},2);
        accSteC(:,iT) = std(accDataC{iT},0,2)./sqrt(nSubjects);
        rtMeanC(:,iT) = mean(rtDataC{iT},2);
        rtSteC(:,iT) = std(rtDataC{iT},0,2)./sqrt(nSubjects);
    end
end

% % all contrasts separately
% for iT = 1:nT % early/late
%     accMeanAll{iT} = squeeze(mean(accDataAll{iT},2));
%     accSteAll{iT} = squeeze(std(accDataAll{iT},0,2)./sqrt(nSubjects));
%     rtMeanAll{iT} = squeeze(mean(rtDataAll{iT},2));
%     rtSteAll{iT} = squeeze(std(rtDataAll{iT},0,2)./sqrt(nSubjects));
% end

p = expt.p;
nCV = numel(p.cueValidity);
tc = p.targetContrasts(contrastIdx)*100;

%% Plot figs
intervalNames = {'T1','T2','T3'};
cueNames = {'valid','invalid','neutral'};
% accLims = [0.4 1];
accLims = [0 3];
rtLims = [0.3 1.4]; 
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
    if iRI==numel(p.respInterval)
        legend(p1, cueNames{idx},'location','SE')
    end
    title(intervalNames{iRI})
    xlim(xlims)
    ylim(accLims)
    rd_supertitle(sprintf('%s run %d, contrast = %d%%, N=%d', exptName, run, tc, nSubjects));
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
    if iRI==numel(p.respInterval)
        legend(cueNames{idx},'location','best')
    end
    title(intervalNames{iRI})
    xlim(xlims)
%     ylim(rtLims)
    ylim([0 1.3])
    box off
    rd_supertitle(sprintf('%s run %d, contrast = %d%%, N=%d', exptName, run, tc, nSubjects));
    rd_raiseAxis(gca);
end

% fig(3) = figure;
% for iRI = 1:numel(p.respInterval)
%     subplot(1,numel(p.respInterval),iRI);
%     hold on
%     
%     plot([0 nCV+1], [0.5 0.5], '--k');
%     b1 = bar(1:nCV, accMean(idx,iRI),'FaceColor',[.5 .5 .5]);
%     p1 = errorbar(1:nCV, accMean(idx,iRI)', accSte(idx,iRI)','k','LineStyle','none');
%     
%     set(gca,'XTick',1:nCV)
% %     set(gca,'XTickLabel', p.cueValidity(idx))
%     set(gca,'XTickLabel', cueNames(idx))
%     xlabel('cue validity')
%     ylabel('acc')
%     title(intervalNames{iRI})
%     ylim([0.3 max(accMean(:))*1.1])
%     ylim([.3 .9])
%     box off
%     rd_supertitle(sprintf('%s run %d, contrast = %d%%, N=%d', exptName, run, tc, nSubjects));
%     rd_raiseAxis(gca);
% end
% 
% fig(4) = figure;
% for iRI = 1:numel(p.respInterval)
%     subplot(1,numel(p.respInterval),iRI);
%     hold on
%     
%     b1 = bar(1:nCV, rtMean(idx,iRI),'FaceColor',[.5 .5 .5]);
%     p1 = errorbar(1:nCV, rtMean(idx,iRI)', rtSte(idx,iRI)','k','LineStyle','none');
%     
%     set(gca,'XTick',1:nCV)
% %     set(gca,'XTickLabel', p.cueValidity(idx))
%     set(gca,'XTickLabel', cueNames(idx))
%     xlabel('cue validity')
%     ylabel('rt')
%     title(intervalNames{iRI})
%     ylim([min(rtMean(:))*0.9 max(rtMean(:))*1.1])
%     box off
%     rd_supertitle(sprintf('%s run %d, contrast = %d%%, N=%d', exptName, run, tc, nSubjects));
%     rd_raiseAxis(gca);
% end

% average across contrasts
fig(5) = figure;
for iRI = 1:numel(p.respInterval)
    subplot(1,numel(p.respInterval),iRI);
    hold on
    
    plot([0 nCV+1], [0.5 0.5], '--k');
    b1 = bar(1:nCV, accMeanC(idx,iRI),'FaceColor',[.5 .5 .5]);
    p1 = errorbar(1:nCV, accMeanC(idx,iRI)', accSteC(idx,iRI)','k','LineStyle','none');
    
    set(gca,'XTick',1:nCV)
%     set(gca,'XTickLabel', p.cueValidity(idx))
    set(gca,'XTickLabel', cueNames(idx))
    xlabel('cue validity')
    ylabel('acc')
    title(intervalNames{iRI})
    ylim([0.3 max(accMeanC(:))*1.1])
%     ylim([.3 .9])
    box off
    rd_supertitle(sprintf('%s run %d, ave across contrasts, N=%d', exptName, run, nSubjects));
    rd_raiseAxis(gca);
end

fig(6) = figure;
for iRI = 1:numel(p.respInterval)
    subplot(1,numel(p.respInterval),iRI);
    hold on
    
    b1 = bar(1:nCV, rtMeanC(idx,iRI),'FaceColor',[.5 .5 .5]);
    p1 = errorbar(1:nCV, rtMeanC(idx,iRI)', rtSteC(idx,iRI)','k','LineStyle','none');
    
    set(gca,'XTick',1:nCV)
%     set(gca,'XTickLabel', p.cueValidity(idx))
    set(gca,'XTickLabel', cueNames(idx))
    xlabel('cue validity')
    ylabel('rt')
    title(intervalNames{iRI})
    ylim([min(rtMeanC(:))*0.9 max(rtMeanC(:))*1.1])
    box off
    rd_supertitle(sprintf('%s run %d, ave across contrasts, N=%d', exptName, run, nSubjects));
    rd_raiseAxis(gca);
end

% % all contrasts separately
% contrastLims = [p.targetContrasts(1)-0.05 p.targetContrasts(end)+0.05];
% fig(7) = figure;
% for iRI = 1:numel(p.respInterval)
%     subplot(1,numel(p.respInterval),iRI)
%     hold on
%     plot(contrastLims, [0.5 0.5], '--k');
%     
%     if numel(p.targetContrasts)>1
%         p1 = errorbar(repmat(p.targetContrasts',1,numel(p.cueValidity)),...
%             accMeanAll{iRI}', accSteAll{iRI}', '.', 'MarkerSize', 20);
%     else
%         for i = 1:length(accMean{iRI})
%             p1(i) = errorbar(p.targetContrasts,...
%                 accMeanAll{iRI}(i), accSteAll{iRI}(i), '.', 'MarkerSize', 20);
%             set(p1(i),'color', colors(i,:))
%         end
%     end
%     xlabel('contrast')
%     ylabel('acc')
%     legend(p1, cueNames,'location','best')
%     title(intervalNames{iRI})
%     xlim(contrastLims)
%     ylim(accLims)
%     rd_supertitle(sprintf('%s run %d, contrast = %d%%, N=%d', exptName, run, tc, nSubjects));
%     rd_raiseAxis(gca);
%     rd_supertitle(axTitle);
% end
% 
% fig(8) = figure;
% rtLims = [0 2];
% for iRI = 1:numel(p.respInterval)
%     subplot(1,numel(p.respInterval),iRI)
%     hold on
%     
%     if numel(p.targetContrasts)>1
%         p1 = errorbar(repmat(p.targetContrasts',1,numel(p.cueValidity)),...
%             rtMeanAll{iRI}', rtSteAll{iRI}', '.', 'MarkerSize', 20);
%     else
%         for i = 1:length(accMean{iRI})
%             p1(i) = errorbar(p.targetContrasts,...
%                 rtMeanAll{iRI}(i), rtSteAll{iRI}(i), '.', 'MarkerSize', 20);
%             set(p1(i),'color', colors(i,:))
%         end
%     end
%     xlabel('contrast')
%     ylabel('RT')
%     legend(p1, cueNames,'location','best')
%     title(intervalNames{iRI})
%     xlim(contrastLims)
%     ylim(rtLims)
%     rd_supertitle(sprintf('%s run %d, contrast = %d%%, N=%d', exptName, run, tc, nSubjects));
%     rd_raiseAxis(gca);
%     rd_supertitle(axTitle);
% end

% invalid breakdown - indiv
fig(9) = figure;
for iRI = 1:numel(p.respInterval)
    subplot(1,numel(p.respInterval),iRI);
    hold on
    plot(xlims, [0.5 0.5], '--k');
    
    p1 = bar(repmat((1:nSubjects)',1,nCV+1),...
        accDataIB{iRI}');
%     colormap(colors(idx,:))
%     colormap(flag(nCV));

    set(gca,'XTick',1:nSubjects)
    set(gca,'XTickLabel', subjectInits)
    xlabel('subject')
    ylabel('acc')
    if iRI==numel(p.respInterval)
        legend(p1, ibNames,'location','SE')
    end
    title(intervalNames{iRI})
    xlim(xlims)
    ylim(accLims)
    rd_supertitle(sprintf('%s run %d, N=%d', exptName, run, nSubjects));
    rd_raiseAxis(gca);
end

fig(10) = figure;
for iRI = 1:numel(p.respInterval)
    subplot(1,numel(p.respInterval),iRI);
    
    p1 = bar(repmat((1:nSubjects)',1,nCV+1),...
        rtDataIB{iRI}');
%     colormap(colors(idx,:))
%     colormap(flag(nCV));
    
    set(gca,'XTick',1:nSubjects)
    set(gca,'XTickLabel', subjectInits)
    xlabel('subject')
    ylabel('rt')
    if iRI==numel(p.respInterval)
        legend(ibNames,'location','best')
    end
    title(intervalNames{iRI})
    xlim(xlims)
%     ylim(rtLims)
    ylim([0 1.3])
    box off
    rd_supertitle(sprintf('%s run %d, N=%d', exptName, run, nSubjects));
    rd_raiseAxis(gca);
end

% invalid breakdown - group
fig(11) = figure;
for iRI = 1:numel(p.respInterval)
    subplot(1,numel(p.respInterval),iRI);
    hold on
    
    plot([0 nCV+2], [0.5 0.5], '--k');
    b1 = bar(1:nCV+1, accMeanIB(:,iRI),'FaceColor',[.5 .5 .5]);
    p1 = errorbar(1:nCV+1, accMeanIB(:,iRI)', accSteIB(:,iRI)','k','LineStyle','none');
    
    set(gca,'XTick',1:nCV+1)
%     set(gca,'XTickLabel', p.cueValidity(idx))
    set(gca,'XTickLabel', ibNames)
    xlabel('cue validity')
    ylabel('acc')
    title(intervalNames{iRI})
    xlim([0 nCV+2])
    ylim([0.3 max(accMeanIB(:))*1.1])
    ylim([.3 .9])
    box off
    rd_supertitle(sprintf('%s run %d, N=%d', exptName, run, nSubjects));
    rd_raiseAxis(gca);
end

fig(12) = figure;
for iRI = 1:numel(p.respInterval)
    subplot(1,numel(p.respInterval),iRI);
    hold on
    
    b1 = bar(1:nCV+1, rtMeanIB(:,iRI),'FaceColor',[.5 .5 .5]);
    p1 = errorbar(1:nCV+1, rtMeanIB(:,iRI)', rtSteIB(:,iRI)','k','LineStyle','none');
    
    set(gca,'XTick',1:nCV+1)
%     set(gca,'XTickLabel', p.cueValidity(idx))
    set(gca,'XTickLabel', ibNames)
    xlabel('cue validity')
    ylabel('rt')
    title(intervalNames{iRI})
    xlim([0 nCV+2])
    ylim([min(rtMeanIB(:))*0.9 max(rtMeanIB(:))*1.1])
    ylim([.2 .8])
    box off
    rd_supertitle(sprintf('%s run %d, N=%d', exptName, run, nSubjects));
    rd_raiseAxis(gca);
end

% 5 bars - group
fig(13) = figure;
for iRI = 1:numel(p.respInterval)
    subplot(1,numel(p.respInterval),iRI);
    hold on
    
    plot([0 nCV+3], [0.5 0.5], '--k');
    b1 = bar(1:nCV, accMeanC(idx,iRI),'FaceColor',[.5 .5 .5]);
    b2 = bar(nCV+1:nCV+2, accMeanIB(3:4,iRI),'FaceColor',[.8 .8 .8]);
    p1 = errorbar(1:nCV, accMeanC(idx,iRI)', accSteC(idx,iRI)','k','LineStyle','none');
    p2 = errorbar(nCV+1:nCV+2, accMeanIB(3:4,iRI)', accSteIB(3:4,iRI)','k','LineStyle','none');
    
    set(gca,'XTick',1:nCV+2)
%     set(gca,'XTickLabel', p.cueValidity(idx))
    set(gca,'XTickLabel', [cueNames(idx) ibNames(3:4)])
    xlabel('cue validity')
    ylabel('acc')
    title(intervalNames{iRI})
    xlim([0 nCV+3])
    ylim([.5 .9])
    box off
    rd_supertitle(sprintf('%s run %d, N=%d', exptName, run, nSubjects));
    rd_raiseAxis(gca);
end

fig(14) = figure;
for iRI = 1:numel(p.respInterval)
    subplot(1,numel(p.respInterval),iRI);
    hold on
    
    b1 = bar(1:nCV, rtMeanC(idx,iRI),'FaceColor',[.5 .5 .5]);
    b2 = bar(nCV+1:nCV+2, rtMeanIB(3:4,iRI),'FaceColor',[.8 .8 .8]);
    p1 = errorbar(1:nCV, rtMeanC(idx,iRI)', rtSteC(idx,iRI)','k','LineStyle','none');
    p2 = errorbar(nCV+1:nCV+2, rtMeanIB(3:4,iRI)', rtSteIB(3:4,iRI)','k','LineStyle','none');
    
    set(gca,'XTick',1:nCV+2)
%     set(gca,'XTickLabel', p.cueValidity(idx))
    set(gca,'XTickLabel', [cueNames(idx) ibNames(3:4)])
    xlabel('cue validity')
    ylabel('rt')
    title(intervalNames{iRI})
    xlim([0 nCV+3])
    ylim([0 .7])
    box off
    rd_supertitle(sprintf('%s run %d, N=%d', exptName, run, nSubjects));
    rd_raiseAxis(gca);
end

%% Set figure properties
for iF = 1:numel(fig)
    % set font size of titles, axis labels, and legends
    set(findall(fig(iF),'type','text'),'FontSize',14)
end

%% Save figs
if saveFigs
    figNames = {'indivAcc','indivRT','groupAcc','groupRT'};
    figPrefix = sprintf('groupData_N%d_%s_run%02d_TemporalAttention3Targets_T1T2all_contrast%d%s', nSubjects, exptName, run, tc, normStr);
    rd_saveAllFigs(fig([1 2 5 6]), figNames, figPrefix)

%     figNames = {'groupAccContrastAve','groupRTContrastAve','groupAccAllContrasts','groupRTAllContrasts'};
%     figPrefix = sprintf('groupData_N%d_%s_run%02d_TemporalAttention3Targets_T1T2all%s', nSubjects, exptName, run, normStr);
%     rd_saveAllFigs(fig(5:end), figNames, figPrefix)
end

%% Reshape data for output to R
acc1 = reshape(accData{1}', nSubjects*3, 1);
acc2 = reshape(accData{2}', nSubjects*3, 1);
acc3 = reshape(accData{3}', nSubjects*3, 1);
rt1 = reshape(rtData{1}', nSubjects*3, 1);
rt2 = reshape(rtData{2}', nSubjects*3, 1);
rt3 = reshape(rtData{3}', nSubjects*3, 1);

acc_all = [acc1; acc2; acc3];
rt_all = [rt1; rt2; rt3];

%% Collapsing across T1/T2/T3
accDataCT = (accDataC{1} + accDataC{2} + accDataC{3})/3;
rtDataCT = (rtDataC{1} + rtDataC{2} + rtDataC{3})/3;

accMeanCT = mean(accDataCT,2);
rtMeanCT = mean(rtDataCT,2);

%% Quick stats
for iT = 1:nT
    fprintf('T%d\n',iT)
    vals = accDataC{iT};
    [hvi pvi cvi svi] = ttest(vals(1,:),vals(2,:));
    [hvn pvn cvn svn] = ttest(vals(1,:),vals(3,:));
    [hni pni cni sni] = ttest(vals(2,:),vals(3,:));
    fprintf('valid vs. invalid, t(%d) = %1.3f, p = %1.4f\n', svi.df, svi.tstat, pvi)
    fprintf('valid vs. neutral, t(%d) = %1.3f, p = %1.4f\n', svn.df, svn.tstat, pvn)
    fprintf('neutral vs. invalid, t(%d) = %1.3f, p = %1.4f\n\n', sni.df, sni.tstat, pni)
end

fprintf('Collapsing across T1, T2, T3\n')
vals = accDataCT;
[hvi pvi cvi svi] = ttest(vals(1,:),vals(2,:));
[hvn pvn cvn svn] = ttest(vals(1,:),vals(3,:));
[hni pni cni sni] = ttest(vals(2,:),vals(3,:));
fprintf('valid vs. invalid, t(%d) = %1.3f, p = %1.4f\n', svi.df, svi.tstat, pvi)
fprintf('valid vs. neutral, t(%d) = %1.3f, p = %1.4f\n', svn.df, svn.tstat, pvn)
fprintf('neutral vs. invalid, t(%d) = %1.3f, p = %1.4f\n\n', sni.df, sni.tstat, pni)

for iT = 1:nT
    fprintf('T%d\n',iT)
    vals = accDataIB{iT};
    [hvi1 pvi1 cvi1 svi1] = ttest(vals(1,:),vals(3,:)); % VI1
    [hvi2 pvi2 cvi2 svi2] = ttest(vals(1,:),vals(4,:)); % VI2
    [hvn pvn cvn svn] = ttest(vals(1,:),vals(2,:)); % VN
    [hni1 pni1 cni1 sni1] = ttest(vals(2,:),vals(3,:)); % NI1
    [hni2 pni2 cni2 sni2] = ttest(vals(2,:),vals(4,:)); % NI2
    fprintf('valid vs. invalid1, t(%d) = %1.3f, p = %1.4f\n', svi1.df, svi1.tstat, pvi1)
    fprintf('valid vs. invalid2, t(%d) = %1.3f, p = %1.4f\n', svi2.df, svi2.tstat, pvi2)
    fprintf('valid vs. neutral, t(%d) = %1.3f, p = %1.4f\n', svn.df, svn.tstat, pvn)
    fprintf('neutral vs. invalid1, t(%d) = %1.3f, p = %1.4f\n', sni1.df, sni1.tstat, pni1)
    fprintf('neutral vs. invalid2, t(%d) = %1.3f, p = %1.4f\n\n', sni2.df, sni2.tstat, pni2)
end

%% Ranomization tests
if doRandomizationTests
    % load empirical null distribution
%     R = load('data/E5_randomizationTest_workspace_run01_N10_20160804.mat');
    R = load('data/E5_randomizationTest_workspace_run01_IB_N10_20160822.mat');
    
    % calculate observed pairwise differences
    accDataCPairwise(1,:) = accMeanC(1,:) - accMeanC(2,:); % VI
    accDataCPairwise(2,:) = accMeanC(1,:) - accMeanC(3,:); % VN
    accDataCPairwise(3,:) = accMeanC(3,:) - accMeanC(2,:); % NI
    
    accDataCTPairwise(1) = accMeanCT(1) - accMeanCT(2);
    accDataCTPairwise(2) = accMeanCT(1) - accMeanCT(3);
    accDataCTPairwise(3) = accMeanCT(3) - accMeanCT(2);
    
    % calculate observed pairwise differences
    rtDataCPairwise(1,:) = rtMeanC(1,:) - rtMeanC(2,:); % VI
    rtDataCPairwise(2,:) = rtMeanC(1,:) - rtMeanC(3,:); % VN
    rtDataCPairwise(3,:) = rtMeanC(3,:) - rtMeanC(2,:); % NI
    
    rtDataCTPairwise(1) = rtMeanCT(1) - rtMeanCT(2);
    rtDataCTPairwise(2) = rtMeanCT(1) - rtMeanCT(3);
    rtDataCTPairwise(3) = rtMeanCT(3) - rtMeanCT(2);
    
    fprintf('RANDOMIZATION TESTS, across contrasts\n')
    fprintf('Accuracy\n')
    for iT = 1:nT
        fprintf('T%d\n',iT)
        for iVC = 1:3 % validity comparison
            pC(iVC,iT) = (nnz(R.accDataCPairwise(iVC,iT,:) > accDataCPairwise(iVC,iT)) + ...
                nnz(R.accDataCPairwise(iVC,iT,:) < -accDataCPairwise(iVC,iT)))/R.nSamples;
        end
        
        fprintf('valid vs. invalid, p = %1.3f\n', pC(1,iT))
        fprintf('valid vs. neutral, p = %1.3f\n', pC(2,iT))
        fprintf('neutral vs. invalid, p = %1.3f\n\n', pC(3,iT))
    end
    
    fprintf('Collapsing across targets\n')
    for iVC = 1:3 % validity comparison
        pCT(iVC) = (nnz(R.accDataCTPairwise(iVC,:) > accDataCTPairwise(iVC)) + ...
            nnz(R.accDataCTPairwise(iVC,:) < -accDataCTPairwise(iVC)))/R.nSamples;
    end
    
    fprintf('valid vs. invalid, p = %1.3f\n', pCT(1))
    fprintf('valid vs. neutral, p = %1.3f\n', pCT(2))
    fprintf('neutral vs. invalid, p = %1.3f\n\n', pCT(3))
    
    fprintf('RT\n')
    for iT = 1:nT
        fprintf('T%d\n',iT)
        for iVC = 1:3 % validity comparison
            pC(iVC,iT) = (nnz(R.rtDataCPairwise(iVC,iT,:) < rtDataCPairwise(iVC,iT)) + ...
                nnz(R.rtDataCPairwise(iVC,iT,:) > -rtDataCPairwise(iVC,iT)))/R.nSamples;
        end
        
        fprintf('valid vs. invalid, p = %1.3f\n', pC(1,iT))
        fprintf('valid vs. neutral, p = %1.3f\n', pC(2,iT))
        fprintf('neutral vs. invalid, p = %1.3f\n\n', pC(3,iT))
    end
    
    fprintf('Collapsing across targets\n')
    for iVC = 1:3 % validity comparison
        pCT(iVC) = (nnz(R.rtDataCTPairwise(iVC,:) < rtDataCTPairwise(iVC)) + ...
            nnz(R.rtDataCTPairwise(iVC,:) > -rtDataCTPairwise(iVC)))/R.nSamples;
    end
    
    fprintf('valid vs. invalid, p = %1.3f\n', pCT(1))
    fprintf('valid vs. neutral, p = %1.3f\n', pCT(2))
    fprintf('neutral vs. invalid, p = %1.3f\n\n', pCT(3))
    
    % invalid breakdown
    % calculate observed pairwise differences
    accDataIBPairwise(1,:) = accMeanIB(1,:) - accMeanIB(3,:); % VI1
    accDataIBPairwise(2,:) = accMeanIB(1,:) - accMeanIB(4,:); % VI2
    accDataIBPairwise(3,:) = accMeanIB(1,:) - accMeanIB(2,:); % VN
    accDataIBPairwise(4,:) = accMeanIB(2,:) - accMeanIB(3,:); % NI1
    accDataIBPairwise(5,:) = accMeanIB(2,:) - accMeanIB(4,:); % NI2
    
    rtDataIBPairwise(1,:) = rtMeanIB(1,:) - rtMeanIB(3,:); % VI1
    rtDataIBPairwise(2,:) = rtMeanIB(1,:) - rtMeanIB(4,:); % VI2
    rtDataIBPairwise(3,:) = rtMeanIB(1,:) - rtMeanIB(2,:); % VN
    rtDataIBPairwise(4,:) = rtMeanIB(2,:) - rtMeanIB(3,:); % NI1
    rtDataIBPairwise(5,:) = rtMeanIB(2,:) - rtMeanIB(4,:); % NI2
    
    fprintf('RANDOMIZATION TESTS, invalid breakdown\n')
    fprintf('Accuracy\n')
    for iT = 1:nT
        fprintf('T%d\n',iT)
        for iVC = 1:5 % validity comparison
            pIB(iVC,iT) = (nnz(R.accDataIBPairwise(iVC,iT,:) > accDataIBPairwise(iVC,iT)) + ...
                nnz(R.accDataIBPairwise(iVC,iT,:) < -accDataIBPairwise(iVC,iT)))/R.nSamples;
        end
        
        fprintf('valid vs. invalid1, p = %1.3f\n', pIB(1,iT))
        fprintf('valid vs. invalid2, p = %1.3f\n', pIB(2,iT))
        fprintf('valid vs. neutral, p = %1.3f\n', pIB(3,iT))
        fprintf('neutral vs. invalid1, p = %1.3f\n', pIB(4,iT))
        fprintf('neutral vs. invalid2, p = %1.3f\n\n', pIB(5,iT))
    end
    
    fprintf('RT\n')
    for iT = 1:nT
        fprintf('T%d\n',iT)
        for iVC = 1:5 % validity comparison
            pIB(iVC,iT) = (nnz(R.rtDataIBPairwise(iVC,iT,:) < rtDataIBPairwise(iVC,iT)) + ...
                nnz(R.rtDataIBPairwise(iVC,iT,:) > -rtDataIBPairwise(iVC,iT)))/R.nSamples;
        end
        
        fprintf('valid vs. invalid1, p = %1.3f\n', pIB(1,iT))
        fprintf('valid vs. invalid2, p = %1.3f\n', pIB(2,iT))
        fprintf('valid vs. neutral, p = %1.3f\n', pIB(3,iT))
        fprintf('neutral vs. invalid1, p = %1.3f\n', pIB(4,iT))
        fprintf('neutral vs. invalid2, p = %1.3f\n\n', pIB(5,iT))
    end
end

%% Effect size
% calculate observed pairwise differences
for iT = 1:nT
    accDataCP(1,:,iT) = accDataC{iT}(1,:) - accDataC{iT}(2,:); % VI
    accDataCP(2,:,iT) = accDataC{iT}(1,:) - accDataC{iT}(3,:); % VN
    accDataCP(3,:,iT) = accDataC{iT}(3,:) - accDataC{iT}(2,:); % NI
end

dP = mean(accDataCP,2)./std(accDataCP,0,2);

for iT = 1:nT
    rtDataCP(1,:,iT) = rtDataC{iT}(1,:) - rtDataC{iT}(2,:); % VI
    rtDataCP(2,:,iT) = rtDataC{iT}(1,:) - rtDataC{iT}(3,:); % VN
    rtDataCP(3,:,iT) = rtDataC{iT}(3,:) - rtDataC{iT}(2,:); % NI
end

dP = mean(rtDataCP,2)./std(rtDataCP,0,2);

% collapsing across target and contrast
accDataCTP(1,:) = accDataCT(1,:) - accDataCT(2,:); % VI
accDataCTP(2,:) = accDataCT(1,:) - accDataCT(3,:); % VN
accDataCTP(3,:) = accDataCT(3,:) - accDataCT(2,:); % NI
 
dCTP = mean(accDataCTP,2)./std(accDataCTP,0,2);

rtDataCTP(1,:) = rtDataCT(1,:) - rtDataCT(2,:); % VI
rtDataCTP(2,:) = rtDataCT(1,:) - rtDataCT(3,:); % VN
rtDataCTP(3,:) = rtDataCT(3,:) - rtDataCT(2,:); % NI
 
dCTP = mean(rtDataCTP,2)./std(rtDataCTP,0,2);

% invalid breakdown
for iT = 1:nT
    accDataIBP(1,:,iT) = accDataIB{iT}(1,:) - accDataIB{iT}(3,:); % VI1
    accDataIBP(2,:,iT) = accDataIB{iT}(1,:) - accDataIB{iT}(4,:); % VI2
    accDataIBP(3,:,iT) = accDataIB{iT}(1,:) - accDataIB{iT}(2,:); % VN
    accDataIBP(4,:,iT) = accDataIB{iT}(2,:) - accDataIB{iT}(3,:); % NI1
    accDataIBP(5,:,iT) = accDataIB{iT}(2,:) - accDataIB{iT}(4,:); % NI2
end

dP = mean(accDataIBP,2)./std(accDataIBP,0,2);

[hh pp] = ttest(accDataIBP,[],[],[],2);

% R: pwr.t.test(d = 1.2335, sig.level = .05, power = .8, type = "paired")
% http://www.statmethods.net/stats/power.html
% http://powerandsamplesize.com/Calculators/Test-1-Mean/1-Sample-Equality




