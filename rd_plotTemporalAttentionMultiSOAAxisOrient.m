% rd_plotTemporalAttentionMultiSOAAxisOrient.m

%% load data
sosa = load('data/E2_SOA_cbD6_run98_N4_SOSA_workspace_20160128.mat');
soda = load('data/E2_SOA_cbD6_run98_N4_SODA_workspace_20160128.mat');
dosa = load('data/E2_SOA_cbD6_run98_N4_DOSA_workspace_20160128.mat');
doda = load('data/E2_SOA_cbD6_run98_N4_DODA_workspace_20160128.mat');

%% get data
m = 'accDataCueAve';
% data for each subject - subject is 3rd dimension
for iT = 1:2
    valsData{iT}(:,1,:) = sosa.(m){iT};
    valsData{iT}(:,2,:) = soda.(m){iT};
    valsData{iT}(:,3,:) = dosa.(m){iT};
    valsData{iT}(:,4,:) = doda.(m){iT};
end
% mean across subjects
for iT = 1:2
    valsMean{iT} = mean(valsData{1},3);
end

%% plot
aoNames = {'sosa','soda','dosa','doda'};
t1t2soa = sosa.t1t2soa;
soaLims = sosa.soaLims;
nSubjects = sosa.nSubjects;
subjectInits = sosa.subjectInits;
accLims = [.4 1];
colors = [0 0 0; 2 132 130; 0 0 0; 2 132 130]/255;

for iSubject = 1:nSubjects
    figure
    for iT = 1:2
        subplot(1,2,iT)
        hold on
        plot(soaLims, [0.5 0.5], '--k');
        set(gca, 'ColorOrder', colors);
        p1 = plot(t1t2soa, valsData{iT}(:,:,iSubject),'-','LineWidth',2);
        
        h = plot(t1t2soa, valsData{iT}(:,1:2,iSubject),'o','LineWidth',1,'MarkerSize',8);
        for iH = 1:numel(h);
            set(h(iH), 'MarkerFaceColor', get(h(iH), 'Color'));
        end
        h = plot(t1t2soa, valsData{iT}(:,3:4,iSubject),'s','LineWidth',1,'MarkerSize',8);
        for iH = 1:numel(h);
            set(h(iH), 'MarkerFaceColor', get(h(iH), 'Color'));
        end
        
        if iT==2
            legend(p1, aoNames)
        end
        xlabel('soa')
        ylabel(m)
        xlim(soaLims)
        ylim(accLims)
    end
    rd_supertitle(sprintf('%s',subjectInits{iSubject}))
end

figure
for iT = 1:2
    subplot(1,2,iT)
    hold on
    plot(soaLims, [0.5 0.5], '--k');
    set(gca, 'ColorOrder', colors);
    p1 = plot(t1t2soa, valsMean{iT},'-','LineWidth',2);
    
    h = plot(t1t2soa, valsMean{iT}(:,1:2),'o','LineWidth',1,'MarkerSize',8);
    for iH = 1:numel(h);
        set(h(iH), 'MarkerFaceColor', get(h(iH), 'Color'));
    end
%     set(h(2), 'MarkerFaceColor', 'w');
    h = plot(t1t2soa, valsMean{iT}(:,3:4),'s','LineWidth',1,'MarkerSize',8);
    for iH = 1:numel(h);
        set(h(iH), 'MarkerFaceColor', get(h(iH), 'Color'));
    end
%     set(h(2), 'MarkerFaceColor', 'w');
    
    
    if iT==2
        legend(p1, aoNames)
    end
    xlabel('soa')
    ylabel(m)
    xlim(soaLims)
    ylim(accLims)
end
rd_supertitle(sprintf('N=%d',nSubjects))
