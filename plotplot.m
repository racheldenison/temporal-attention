% plotplot.m

% colors = {'b','r',[.5 .5 .5]};
% fcolors = {'b','r',[.5 .5 .5]};
colors = {[0 143 133]/255,[212 36 113]/255,[.5 .5 .5]};
fcolors = colors;
dpLims = [-0.5 3];

fig(1) = figure;
for iRI = 1:numel(p.respInterval)
    subplot(1,numel(p.respInterval),iRI)
    hold on
    plot(soaLims, [0 0], '--k');
    
    for iCV = [3 2 1]
        
        p1 = plot(t1t2soa,...
            dpMean{iRI}(iCV,:), '-','LineWidth', 1.5, ...
            'Color', colors{iCV});
        e1 = errorbar(t1t2soa,...
            dpMean{iRI}(iCV,:), dpSte{iRI}(iCV,:), 'o', ...
            'LineWidth', 1, 'Color', colors{iCV}, ...
            'MarkerSize',8,'MarkerFaceColor',fcolors{iCV});
    end
    
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

% colors = {[102,45,145]/255, [241,90,41]/255};
colors = {[0,102,204]/255, [228,26,28]/255};
fcolors = colors;
fcolors = {'w','w','w'};
fig(2) = figure;
hold on
plot(soaLims, [0 0], '--k')
p1 = [];
p1(1) = plot(t1t2soa, mean(accDataCueEff{1},2),'-','LineWidth',2.1, 'Color', colors{1});
p1(2) = plot(t1t2soa, mean(accDataCueEff{2},2),'-','LineWidth',2.1, 'Color', colors{2});

e1(1) = errorbar(t1t2soa, mean(accDataCueEff{1},2), ...
    std(accDataCueEff{1},0,2)./sqrt(nSubjects),...
    'o','MarkerSize', 8,'LineWidth', 1.1, 'Color', colors{1}, 'MarkerFaceColor',fcolors{1});
e1(2) = errorbar(t1t2soa, mean(accDataCueEff{2},2), ...
    std(accDataCueEff{2},0,2)./sqrt(nSubjects),...
    'o','MarkerSize', 8,'LineWidth', 1.1, 'Color', colors{2}, 'MarkerFaceColor',fcolors{2});

legend(p1, intervalNames)
xlabel('soa')
ylabel('cuing effect (accuracy valid-invalid)')
xlim(soaLims)
% ylim(accCueEffLims)
ylim([-.4 .8])

colors = {[0 143 133]/255,[212 36 113]/255,[.5 .5 .5]};
fcolors = colors;

fig(3) = figure;
for iRI = 1:numel(p.respInterval)
    subplot(3,numel(p.respInterval),[iRI iRI+2])
    hold on
    plot(soaLims, [0 0], '--k');
    
    for iCV = [3 2 1]
        
        p1 = plot(t1t2soa,...
            dpMean{iRI}(iCV,:), '-','LineWidth', 1.5, ...
            'Color', colors{iCV});
        e1 = errorbar(t1t2soa,...
            dpMean{iRI}(iCV,:), dpSte{iRI}(iCV,:), 'o', ...
            'LineWidth', 1, 'Color', colors{iCV}, ...
            'MarkerSize',8,'MarkerFaceColor',fcolors{iCV});
    end
    
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

for iRI = 1:numel(p.respInterval)
    subplot(3,numel(p.respInterval),iRI+4);
    hold on
    plot(soaLims, [0 0], '--k')
    p1 = plot(t1t2soa, mean(accDataCueEff{iRI},2),'-','LineWidth',2.1, 'Color', 'k');
    e1 = errorbar(t1t2soa, mean(accDataCueEff{iRI},2), ...
    std(accDataCueEff{1},0,2)./sqrt(nSubjects),...
        'o','MarkerSize', 8,'LineWidth', 1.1, 'Color', 'k', 'MarkerFaceColor','w');
    xlim(soaLims)
    ylim([-.4 .8])
end
