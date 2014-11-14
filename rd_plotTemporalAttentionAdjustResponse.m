% rd_plotTemporalAttentionAdjustResponse.m

%% setup
subject = 'bl_a1_tc100_soa1000-1250';
run = 9;

smoothHist = 0;

expName = 'E3_adjust';
% dataDir = 'data';
% figDir = 'figures';
dataDir = pathToExpt('data');
figDir = pathToExpt('figures');
dataDir = sprintf('%s/%s/%s', dataDir, expName, subject(1:2));
figDir = sprintf('%s/%s/%s', figDir, expName, subject(1:2));

binSize = 4;
edges = -90:binSize:90;

%% load data
dataFile = dir(sprintf('%s/%s_run%02d*', dataDir, subject, run));
load(sprintf('%s/%s', dataDir, dataFile.name))

%% analyze and plot
errorIdx = strcmp(expt.trials_headers, 'responseError');

% for iEL = 1:2
%     for iV = 1:3
%         count = histc(results.totals.all{iV,iEL}(:,errorIdx),edges);
%         nTrials = size(results.totals.all{iV,iEL},1);
%         rate{iEL}(iV,:) = count./nTrials;
%     end
% end

for iEL = 1:2
    for iV = 1:3
        vals = results.totals.all{iV,iEL}(:,errorIdx);
        nTrials = size(vals,1);
        
        if smoothHist
            counts = nan(binSize*2, (edges(end)+binSize) - (edges(1)-binSize));
            for iCount = 1:binSize*2+1
                count = histc(vals, edges-binSize+iCount-1);
                counts(iCount,(1:binSize:size(counts,2)-binSize)+iCount-1) = count;
            end
            rate0 = nanmean(counts,1)./nTrials;
            rate{iEL}(iV,:) = rate0(binSize+1:end-binSize);
            xEdges = edges(1):edges(end);
        else
            count = histc(results.totals.all{iV,iEL}(:,errorIdx),edges);
            rate{iEL}(iV,:) = count./nTrials;
            xEdges = edges;
        end
    end
end

colors = {'b','g','r'};
figure
for iEL = 1:2
    subplot(1,2,iEL)
    hold on
    for iV = 1:3
        vals = rate{iEL}(iV,:);
        a1 = area(xEdges, vals, 'FaceColor', colors{iV}, ...
            'EdgeColor', colors{iV}, 'LineWidth', 1.5);
        child = get(a1, 'Children');
        set(child,'FaceAlpha',.5);
    end
    xlim([-90 90])
    ylim([0 .5])
    xlabel('error in orientation report')
    ylabel('p(error)')
end
legend('valid','invalid','neutral')
legend('boxoff')