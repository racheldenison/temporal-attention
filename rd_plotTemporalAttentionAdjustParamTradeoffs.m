% rd_plotTemporalAttentionAdjustParamTradeoffs.m
%
% plot tradeoffs between fit parameters, eg. sd and g

%% setup
% see rd_organizeAdjustFitGroupStats.m
load data/adjust_fit_group_stats_mixtureNoBiasMaxPosterior_run09_N12_20150512.mat

nSubjects = numel(subjectIDs);

targetNames = {'T1','T2'};
validityNames = {'valid','invalid','neutral'};

measures = fields(paramsData);

nVal = size(paramsData.sd,1);
nTarg = size(paramsData.sd,2);
nSub = size(paramsData.sd,3);

%% plot tradeoffs (should not really go in this script)
sdT1 = squeeze(paramsData.sd(:,1,:));
sdT2 = squeeze(paramsData.sd(:,2,:));
gT1 = squeeze(paramsData.g(:,1,:));
gT2 = squeeze(paramsData.g(:,2,:));

xlims = [-16 16];
ylims = [-.2 .2];

figure
subplot(1,2,1)
hold on
plot(sdT1(2,:)-sdT1(1,:), gT1(2,:)-gT1(1,:),'.')
xlim(xlims)
ylim(ylims)
vline(0,'k')
plot(xlims,[0 0],'k')
axis square
xlabel('standard deviation (invalid-valid)')
ylabel('guess rate (invalid-valid)')
title('T1')
subplot(1,2,2)
hold on
plot(sdT2(2,:)-sdT2(1,:), gT2(2,:)-gT2(1,:),'.')
xlim(xlims)
ylim(ylims)
vline(0,'k')
plot(xlims,[0 0],'k')
axis square
title('T2')

[rho p] = corr((sdT1(2,:)-sdT1(1,:))', (gT1(2,:)-gT1(1,:))');
[rho p] = corr((sdT2(2,:)-sdT2(1,:))', (gT2(2,:)-gT2(1,:))');