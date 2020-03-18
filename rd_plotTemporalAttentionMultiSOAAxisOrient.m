% rd_plotTemporalAttentionMultiSOAAxisOrient.m

%% load data
% sosa = load('data/E2_SOA_cbD6_run98_N4_SOSA_workspace_20160128.mat');
% soda = load('data/E2_SOA_cbD6_run98_N4_SODA_workspace_20160128.mat');
% dosa = load('data/E2_SOA_cbD6_run98_N4_DOSA_workspace_20160128.mat');
% doda = load('data/E2_SOA_cbD6_run98_N4_DODA_workspace_20160128.mat');

% sosa = load('data/E2_SOA_cbD6_run98_N5_SOSA_workspace_20180731.mat');
% soda = load('data/E2_SOA_cbD6_run98_N5_SODA_workspace_20180731.mat');
% dosa = load('data/E2_SOA_cbD6_run98_N5_DOSA_workspace_20180731.mat');
% doda = load('data/E2_SOA_cbD6_run98_N5_DODA_workspace_20180731.mat');

sosa = load('data/E2_SOA_cbD6_run98_N5_SOSA_workspace_20191009.mat');
soda = load('data/E2_SOA_cbD6_run98_N5_SODA_workspace_20191009.mat');
dosa = load('data/E2_SOA_cbD6_run98_N5_DOSA_workspace_20191009.mat');
doda = load('data/E2_SOA_cbD6_run98_N5_DODA_workspace_20191009.mat');

%% get data
m = 'dpData'; % dpData, dpDataCueEff
% data for each subject - subject is 3rd dimension
for iT = 1:2
    valsData{iT}(:,1,:) = sosa.(m){iT}(2,:,:); %(1,:,:) for dpData, or none for dpDataCueEff;
    valsData{iT}(:,2,:) = soda.(m){iT}(2,:,:); %(1,:,:);
    valsData{iT}(:,3,:) = dosa.(m){iT}(2,:,:); %(1,:,:);
    valsData{iT}(:,4,:) = doda.(m){iT}(2,:,:); %(1,:,:);
end
% mean across subjects
for iT = 1:2
    valsMean{iT} = mean(valsData{iT},3);
end

%% plot
aoNames = {'sosa','soda','dosa','doda'};
t1t2soa = sosa.t1t2soa;
soaLims = sosa.soaLims;
nSubjects = sosa.nSubjects;
subjectInits = sosa.subjectInits;
valLims = [-0.5 3.5]; %[.4 1];
colors = [0 0 0; 2 132 130; 0 0 0; 2 132 130]/255;

for iSubject = 1:nSubjects
    figure
    for iT = 1:2
        subplot(1,2,iT)
        hold on
        plot(soaLims, [0 0], '--k');
        set(gca, 'ColorOrder', colors);
        p1 = plot(t1t2soa, valsData{iT}(:,:,iSubject),'-','LineWidth',2);
        
        h = plot(t1t2soa, valsData{iT}(:,1:2,iSubject),'o','LineWidth',1,'MarkerSize',8);
        for iH = 1:numel(h)
            set(h(iH), 'MarkerFaceColor', get(h(iH), 'Color'));
        end
        h = plot(t1t2soa, valsData{iT}(:,3:4,iSubject),'o','LineWidth',1,'MarkerSize',8);
        for iH = 1:numel(h)
%             set(h(iH), 'MarkerFaceColor', get(h(iH), 'Color'));
            set(h(iH), 'MarkerFaceColor', 'w');
        end
        
        if iT==2
            legend(p1, aoNames)
        end
        xlabel('soa')
        ylabel(m)
        xlim(soaLims)
        ylim(valLims)
    end
    rd_supertitle(sprintf('%s',subjectInits{iSubject}))
end

figure
for iT = 1:2
    subplot(1,2,iT)
    hold on
    plot(soaLims, [0 0], '--k');
    set(gca, 'ColorOrder', colors);
    p1 = plot(t1t2soa, valsMean{iT},'-','LineWidth',2);
    
    h = plot(t1t2soa, valsMean{iT}(:,1:2),'o','LineWidth',1,'MarkerSize',8);
    for iH = 1:numel(h)
        set(h(iH), 'MarkerFaceColor', get(h(iH), 'Color'));
    end
    allh = h;
    h = plot(t1t2soa, valsMean{iT}(:,3:4),'o','LineWidth',1,'MarkerSize',8);
    for iH = 1:numel(h)
%         set(h(iH), 'MarkerFaceColor', get(h(iH), 'Color'));
        set(h(iH), 'MarkerFaceColor', 'w');
    end
	allh = [allh; h];
    
    
    if iT==2
        legend(allh, aoNames)
    end
    xlabel('soa')
    ylabel(m)
    xlim(soaLims)
    ylim(valLims)
end
rd_supertitle(sprintf('N=%d',nSubjects))
