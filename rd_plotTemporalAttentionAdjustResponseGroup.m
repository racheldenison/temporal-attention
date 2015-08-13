% rd_plotTemporalAttentionAdjustResponseGroup.m

%% setup
subjectIDs = {'bl','rd','id','ec','ld','en','sj','ml','ca','jl','ew','jx'};
run = 9;
nSubjects = numel(subjectIDs);

saveFigs = 0;

%% get data
for iSubject = 1:nSubjects
    subjectID = subjectIDs{iSubject};
    [groupData.error(iSubject,:), groupData.pError(iSubject,:)] = ...
        rd_plotTemporalAttentionAdjustResponse(subjectID, run);
end
xEdges = mean(groupData.error,1);

for iSubject = 1:nSubjects
    for iEL = 1:2
        groupData.pError2{iEL}(:,:,iSubject) = groupData.pError{iSubject,iEL};
        groupMean.rate{iEL} = mean(groupData.pError2{iEL},3);
        groupSte.rate{iEL} = std(groupData.pError2{iEL},0,3)./sqrt(nSubjects);
    end
end

%% plot
%% overlaid
plotType = 'shadedError'; % 'area','line','shadedError'
colors = {'b','g','r'};
fig(1) = figure;
for iEL = 1:2
    subplot(1,2,iEL)
    hold on
    for iV = 1:3
        vals = groupMean.rate{iEL}(iV,:);
        valsSte = groupSte.rate{iEL}(iV,:);
        switch plotType
            case 'area'
                a1 = area(xEdges, vals, 'FaceColor', colors{iV}, ...
                    'EdgeColor', colors{iV}, 'LineWidth', 1.5);
                child = get(a1, 'Children');
                set(child,'FaceAlpha',.5);
            case 'line'
                p1 = plot(xEdges, vals, 'Color', colors{iV}, ...
                    'LineWidth', 1);
            case 'shadedError'
                p1 = shadedErrorBar(xEdges, vals, valsSte, colors{iV}, 1);
            otherwise
                error('plotType not recognized')
        end
    end
    xlim([-90 90])
    ylim([0 .08])
    xlabel('error in orientation report')
    ylabel('p(error)')
    
    if iEL == 2
        legend('valid','invalid','neutral')
        legend('boxoff')
    end
    plot_vertical_line(0);
end

rd_supertitle(sprintf('%s ', subjectIDs{:}));
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
        vals = groupMean.rate{iEL}(v,:);
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
rd_supertitle(sprintf('%s ', subjectIDs{:}));
rd_raiseAxis(gca);

%% save figs
if saveFigs
    figNames = {'errorHistOverlay','errorHistSeparate'};
    figPrefix = sprintf('gE3_N%d_run%02d', nSubjects, run);
    rd_saveAllFigs(fig, figNames, figPrefix, [], '-depsc2')
end
