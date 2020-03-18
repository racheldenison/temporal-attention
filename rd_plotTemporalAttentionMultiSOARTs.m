% rd_plotTemporalAttentionMultiSOARTs.m

%% setup
%     dataDir = '/Users/rachel/Documents/NYU/Projects/Temporal_Attention/Code/Expt_Scripts/Behav/data';
    dataDir = '/Local/Users/denison/Google Drive/NYU/Projects/Temporal_Attention/Code/Expt_Scripts/Behav/data';

subjects = {'rd','hl','ho','vp','jp'};
nSubjects = numel(subjects);
nSOA = 10;

%% take a bootstrap from trial data aggregated across subjects
allRT = [];
for iS = 1:nSubjects
    rt = [];
    for iSOA = 1:nSOA
        fileName = sprintf('%s/E2_SOA_cbD6_%s_run98_alldata', dataDir, subjects{iS});
        load(fileName)
        
        expt = E{iSOA};
        rt(:,iSOA) = expt.trials(:,strcmp(expt.trials_headers,'rt'));
    end
    allRT{iS} = rt(:);
end

%% how many trials longer than a certain threshold?
threshs = [1 2];
for iS = 1:nSubjects
    for iTh = 1:numel(threshs)
        pLong(iTh,iS) = nnz(allRT{iS}>threshs(iTh))/numel(allRT{iS});
    end
end

%% plot 
% RT histograms
figure
hold on
for iS = 1:nSubjects
    histogram(allRT{iS},'Normalization','pdf','DisplayStyle','stairs')
end
xlabel('RT (s)')
ylabel('pdf')
legend(subjects)

% percent RT above a certain thresh
figure
bar(pLong*100)
set(gca,'XTickLabel',threshs)
xlabel('RT threshold (s)')
ylabel('% > threshold')
legend(subjects)