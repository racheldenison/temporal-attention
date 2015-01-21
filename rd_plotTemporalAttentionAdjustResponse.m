function [xEdges, rate] = rd_plotTemporalAttentionAdjustResponse(subjectID, run)

%% setup
% subjectID = 'xx';
subject = sprintf('%s_a1_tc100_soa1000-1250', subjectID);
% run = 9;

saveFigs = 1;

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
binCenters = edges(1:end-1) + binSize/2;

%% load data
dataFile = dir(sprintf('%s/%s_run%02d*', dataDir, subject, run));
load(sprintf('%s/%s', dataDir, dataFile(1).name))

%% analyze and plot
errorIdx = strcmp(expt.trials_headers, 'responseError');

for iEL = 1:2
    for iV = 1:3
        vals = results.totals.all{iV,iEL}(:,errorIdx);
        nTrials = size(vals,1);
        
        if smoothHist
            % n.b. counts not recentered to bin centers, and rate not
            % adjusted for bin size -- best not to use at this stage
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
            count(end-1) = count(end-1) + count(end); % merge last bin with previous
            count(end) = [];
%             rate{iEL}(iV,:) = count./nTrials;
            rate{iEL}(iV,:) = count./sum(count*binSize);
            xEdges = binCenters;
        end
    end
end

%% overlaid
plotType = 'area'; % 'area','line'
colors = {'b','g','r'};
fig(1) = figure;
for iEL = 1:2
    subplot(1,2,iEL)
    hold on
    for iV = 1:3
        vals = rate{iEL}(iV,:);
        switch plotType
            case 'area'
                a1 = area(xEdges, vals, 'FaceColor', colors{iV}, ...
                    'EdgeColor', colors{iV}, 'LineWidth', 1.5);
                child = get(a1, 'Children');
                set(child,'FaceAlpha',.5);
            case 'line'
                p1 = plot(xEdges, vals, 'Color', colors{iV}, ...
                    'LineWidth', 1);
            otherwise
                error('plotType not recognized')
        end
    end
    xlim([-90 90])
    ylim([0 .08])
    xlabel('error in orientation report')
    ylabel('p(error)')
end
legend('valid','invalid','neutral')
legend('boxoff')
rd_supertitle(sprintf('%s run %d', subjectID, run));
rd_raiseAxis(gca);

%% separate plots for each validity
fig(2) = figure;
targetNames = {'T1','T2'};
validityNames = {'valid','invalid','neutral'};
validityOrder = [1 3 2];
for iEL = 1:2
    for iV = 1:3
        v = validityOrder(iV);
        subplot(3,2,2*(iV-1)+iEL)
        hold on
        plot_vertical_line(0);
        vals = rate{iEL}(v,:);
        switch plotType
            case 'area'
                a1 = area(xEdges, vals, 'FaceColor', colors{v}, ...
                    'EdgeColor', colors{v}, 'LineWidth', 1.5);
                child = get(a1, 'Children');
                set(child,'FaceAlpha',.5);
            case 'line'
                p1 = plot(xEdges, vals, 'Color', colors{v}, ...
                    'LineWidth', 1);
            otherwise
                error('plotType not recognized')
        end
        xlim([-90 90])
        ylim([0 .08])
        if iEL==2
            legend(a1, validityNames{v})
            legend('boxoff')
        end
        if iV==1
            title(targetNames{iEL})
        end
    end
    xlabel('error in orientation report')
    ylabel('p(error)')
end
rd_supertitle(sprintf('%s run %d', subjectID, run));
rd_raiseAxis(gca);

%% save figs
if saveFigs
    figNames = {'errorHistOverlay','errorHistSeparate'};
    figPrefix = sprintf('%s_run%02d_TemporalAttentionAdjust', subject, run);
    rd_saveAllFigs(fig, figNames, figPrefix, figDir, '-depsc2')
end
    
