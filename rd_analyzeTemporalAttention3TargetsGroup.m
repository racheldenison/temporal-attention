% rd_analyzeTemporalAttention3TargetsGroup.m

exptName = 'cbD15';
subjectInits = {'gb','xw','yz','jg','rd'};
contrastIdx = 1;

run = 1;

normalizeData = 0;

saveFigs = 0;

nSubjects = numel(subjectInits);

% dataDir = pathToExpt('data');
% dataDir = [pathToExpt('data') '/pilot/rd'];

doRandomizationTests = 0; % requries having generated empirical null distributions

%% Get data
for iSubject = 1:nSubjects 
    subjectInit = subjectInits{iSubject};
    
    dataDir = sprintf('%s/E5_T3_cbD15/%s', pathToExpt('data'), subjectInit(1:2));
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
%     [expt results] = rd_analyzeTemporalAttention(expt, 0, 0, 0, 0, T1T2Axis, 0);
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    nT = numel(expt.p.respInterval);

    % read out the accuracy and rt
    for iT = 1:nT % early/late
        accData{iT}(:,iSubject) = results.accMean{iT}(:,contrastIdx);
        rtData{iT}(:,iSubject) = results.rtMean{iT}(:,contrastIdx);
    end
    
    % also gather all the data for all contrasts
    for iContrast = 1:numel(expt.p.targetContrasts)
        for iT = 1:nT % early/late
            accDataAll{iT}(:,iSubject,iContrast) = results.accMean{iT}(:,iContrast);
            rtDataAll{iT}(:,iSubject,iContrast) = results.rtMean{iT}(:,iContrast);
            
            % average across contrasts
            accDataC{iT}(:,iSubject) = mean(results.accMean{iT},2);
            rtDataC{iT}(:,iSubject) = mean(results.rtMean{iT},2);
        end
    end
end

% summary across subjects
for iT = 1:nT % early/late
    accMean(:,iT) = mean(accData{iT},2);
    accSte(:,iT) = std(accData{iT},0,2)./sqrt(nSubjects);
    rtMean(:,iT) = mean(rtData{iT},2);
    rtSte(:,iT) = std(rtData{iT},0,2)./sqrt(nSubjects);
end

% average across contrasts
if normalizeData
    % use the Morey 2008 correction
    a = cat(3, accDataC{1}, accDataC{2});
    b = normalizeDC(shiftdim(a,2));
    c = shiftdim(b,1);
    accDataC{1} = c(:,:,1);
    accDataC{2} = c(:,:,2);
    
    a = cat(3, rtDataC{1}, rtDataC{2});
    b = normalizeDC(shiftdim(a,2));
    c = shiftdim(b,1);
    rtDataC{1} = c(:,:,1);
    rtDataC{2} = c(:,:,2);
    
    [fixed1 fixed2 N] = size(b);
    M = fixed1*fixed2;
    morey = M/(M-1);
    
    for iT = 1:nT
        accMeanC(:,iT) = mean(accDataC{iT},2);
        accSteC(:,iT) = sqrt(morey*var(accDataC{iT},0,2)./(nSubjects));
        rtMeanC(:,iT) = mean(rtDataC{iT},2);
        rtSteC(:,iT) = sqrt(morey*var(rtDataC{iT},0,2)./(nSubjects));
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

% all contrasts separately
for iT = 1:nT % early/late
    accMeanAll{iT} = squeeze(mean(accDataAll{iT},2));
    accSteAll{iT} = squeeze(std(accDataAll{iT},0,2)./sqrt(nSubjects));
    rtMeanAll{iT} = squeeze(mean(rtDataAll{iT},2));
    rtSteAll{iT} = squeeze(std(rtDataAll{iT},0,2)./sqrt(nSubjects));
end

p = expt.p;
nCV = numel(p.cueValidity);
tc = p.targetContrasts(contrastIdx)*100;

%% Plot figs
intervalNames = {'T1','T2','T3'};
cueNames = {'valid','invalid','neutral'};
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
    legend(p1, cueNames,'location','best')
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
    legend(cueNames,'location','best')
    title(intervalNames{iRI})
    xlim(xlims)
%     ylim(rtLims)
    ylim([0.2 4])
    box off
    rd_supertitle(sprintf('%s run %d, contrast = %d%%, N=%d', exptName, run, tc, nSubjects));
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
%     set(gca,'XTickLabel', p.cueValidity(idx))
    set(gca,'XTickLabel', cueNames(idx))
    xlabel('cue validity')
    ylabel('acc')
    title(intervalNames{iRI})
    ylim([0.3 max(accMean(:))*1.1])
    ylim([.3 .9])
    box off
    rd_supertitle(sprintf('%s run %d, contrast = %d%%, N=%d', exptName, run, tc, nSubjects));
    rd_raiseAxis(gca);
end

fig(4) = figure;
for iRI = 1:numel(p.respInterval)
    subplot(1,numel(p.respInterval),iRI);
    hold on
    
    b1 = bar(1:nCV, rtMean(idx,iRI),'FaceColor',[.5 .5 .5]);
    p1 = errorbar(1:nCV, rtMean(idx,iRI)', rtSte(idx,iRI)','k','LineStyle','none');
    
    set(gca,'XTick',1:nCV)
%     set(gca,'XTickLabel', p.cueValidity(idx))
    set(gca,'XTickLabel', cueNames(idx))
    xlabel('cue validity')
    ylabel('rt')
    title(intervalNames{iRI})
    ylim([min(rtMean(:))*0.9 max(rtMean(:))*1.1])
    box off
    rd_supertitle(sprintf('%s run %d, contrast = %d%%, N=%d', exptName, run, tc, nSubjects));
    rd_raiseAxis(gca);
end

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
    ylim([.3 .9])
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

%% Set figure properties
for iF = 1:numel(fig)
    % set font size of titles, axis labels, and legends
    set(findall(fig(iF),'type','text'),'FontSize',14)
end

%% Save figs
if saveFigs
    figNames = {'indivAcc','indivRT','groupAcc','groupRT'};
    figPrefix = sprintf('groupData_N%d_%s_run%02d_TemporalAttention3Targets_T1T2all_contrast%d', nSubjects, exptName, run, tc);
    rd_saveAllFigs(fig(1:4), figNames, figPrefix, [], '-depsc2')

    figNames = {'groupAccContrastAve','groupRTContrastAve','groupAccAllContrasts','groupRTAllContrasts'};
    figPrefix = sprintf('groupData_N%d_%s_run%02d_TemporalAttention3Targets_T1T2all', nSubjects, exptName, run);
    rd_saveAllFigs(fig(5:end), figNames, figPrefix, [], '-depsc2')
end

%% Reshape data for output to R
acc1 = reshape(accData{1}', nSubjects*3, 1);
acc2 = reshape(accData{2}', nSubjects*3, 1);
rt1 = reshape(rtData{1}', nSubjects*3, 1);
rt2 = reshape(rtData{2}', nSubjects*3, 1);

acc_all = [acc1; acc2];
rt_all = [rt1; rt2];

%% Collapsing across T1/T2
accDataCT = (accDataC{1} + accDataC{2})/2;
rtDataCT = (rtDataC{1} + rtDataC{2})/2;

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

fprintf('Collapsing across T1 and T2\n')
vals = accDataCT;
[hvi pvi cvi svi] = ttest(vals(1,:),vals(2,:));
[hvn pvn cvn svn] = ttest(vals(1,:),vals(3,:));
[hni pni cni sni] = ttest(vals(2,:),vals(3,:));
fprintf('valid vs. invalid, t(%d) = %1.3f, p = %1.4f\n', svi.df, svi.tstat, pvi)
fprintf('valid vs. neutral, t(%d) = %1.3f, p = %1.4f\n', svn.df, svn.tstat, pvn)
fprintf('neutral vs. invalid, t(%d) = %1.3f, p = %1.4f\n\n', sni.df, sni.tstat, pni)

%% Ranomization tests
if doRandomizationTests
    % load empirical null distribution
    R = load('data/E0_randomizationTest_workspace_run09_N10_20160107b.mat');
    
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
    
    fprintf('Collapsing across T1 and T2\n')
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
    
    fprintf('Collapsing across T1 and T2\n')
    for iVC = 1:3 % validity comparison
        pCT(iVC) = (nnz(R.rtDataCTPairwise(iVC,:) < rtDataCTPairwise(iVC)) + ...
            nnz(R.rtDataCTPairwise(iVC,:) > -rtDataCTPairwise(iVC)))/R.nSamples;
    end
    
    fprintf('valid vs. invalid, p = %1.3f\n', pCT(1))
    fprintf('valid vs. neutral, p = %1.3f\n', pCT(2))
    fprintf('neutral vs. invalid, p = %1.3f\n\n', pCT(3))
end

%% Effect size
% calculate observed pairwise differences
for iT = 1:nT
    accDataCP(1,:,iT) = accDataC{iT}(1,:) - accDataC{iT}(2,:); % VI
    accDataCP(2,:,iT) = accDataC{iT}(1,:) - accDataC{iT}(3,:); % VN
    accDataCP(3,:,iT) = accDataC{iT}(3,:) - accDataC{iT}(2,:); % NI
end

dP = mean(accDataCP,2)./std(accDataCP,0,2);

% R: pwr.t.test(d = 1.2335, sig.level = .05, power = .8, type = "paired")
% http://www.statmethods.net/stats/power.html
% http://powerandsamplesize.com/Calculators/Test-1-Mean/1-Sample-Equality




