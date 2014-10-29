% rd_plotTemporalAttentionAdjustResponse.m

binSize = 2;
edges = 0:binSize:90;

% load data
load([pathToExpt '/data/pilot/rd/rdTest_run09_TemporalAttentionAdjust_T1T2all_20141029.mat'])

errorIdx = strcmp(expt.trials_headers, 'responseError');

for iEL = 1:2
    for iV = 1:3
        count = histc(results.totals.all{iV,iEL}(:,11),edges);
        nTrials = size(results.totals.all{iV,iEL},1);
        rate{iEL}(iV,:) = count./nTrials;
    end
end

colors = {'b','g','r'};
figure
for iEL = 1:2
    subplot(1,2,iEL)
    hold on
    for iV = 1:3
        vals = rate{iEL}(iV,:);
        a1 = area(edges, vals, 'FaceColor', colors{iV}, ...
            'EdgeColor', colors{iV}, 'LineWidth', 1.5);
        child = get(a1, 'Children');
        set(child,'FaceAlpha',.5);
    end
    xlim([0 90])
    ylim([0 .5])
    xlabel('error in orientation report')
    ylabel('p(error)')
end
legend('valid','invalid','neutral')
legend('boxoff')