% rd_plotTemporalAttentionAdjustResponseGroup.m

%% setup
subjectIDs = {'bl','rd','id'};
nSubjects = numel(subjectIDs);

%% get data
for iSubject = 1:nSubjects
    subjectID = subjectIDs{iSubject};
    [groupData.error(iSubject,:), groupData.pError(iSubject,:)] = ...
        rd_plotTemporalAttentionAdjustResponse(subjectID);
end
xEdges = mean(groupData.error,1);

for iSubject = 1:nSubjects
    for iEL = 1:2
        groupData.pError2{iEL}(:,:,iSubject) = groupData.pError{iSubject,iEL};
        rate{iEL} = mean(groupData.pError2{iEL},3);
    end
end

%% plot
%% overlaid
plotType = 'area'; % 'area','line'
colors = {'b','g','r'};
figure
for iEL = 1:2
    subplot(1,2,iEL)
    hold on
    plot_vertical_line(0);
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
    ylim([0 .18])
    xlabel('error in orientation report')
    ylabel('p(error)')
end
legend('valid','invalid','neutral')
legend('boxoff')

%% separate plots for each validity
figure
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
        ylim([0 .18])
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

