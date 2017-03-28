% rd_plotMainEffectsTemporalAttention.m

%% setup
e0 = load('data/E0_workspace_run09_N10_20160224.mat');
e3 = load('data/E3_workspace_run09_N12_20160224.mat');
e5 = load('data/E5_workspace_run01_N10_20160806.mat');

expNames = {'e0','e3','e5'};
nExp = numel(expNames);

normalizeData = 1;

measureOption = 'rt-g'; % 'acc-sd','rt-g'

%% get data
switch measureOption
    case 'acc-sd'
        groupData.e0 = e0.accDataCT;
        groupData.e3 = squeeze(mean(e3.paramsData.sd,2));
        groupData.e5 = e5.accDataCT;
    case 'rt-g'
        groupData.e0 = e0.rtDataCT;
        groupData.e3 = squeeze(mean(e3.paramsData.g,2));
        groupData.e5 = e5.rtDataCT;
    otherwise
        error('measureOption not recognized')
end

%% summarize data, with normalization option
for iExp = 1:nExp
    expName = expNames{iExp};
    nSubjects = size(groupData.(expName),2);
    groupMean.(expName) = mean(groupData.(expName),2);

    if normalizeData
        % use the Morey 2008 correction
        vals = groupData.(expName);
        valsn = normalizeDC(vals);
        
        [fixed(1) N] = size(valsn);
        M = prod(fixed);
        morey = M/(M-1);
        
        groupSte.(expName) = sqrt(morey*var(valsn,0,2)./(nSubjects));
    else
        groupSte.(expName) = std(groupData.(expName),0,2)/sqrt(nSubjects);
    end
end

%% plot
switch measureOption
    case 'acc-sd'
        ylims.e0 = [.5 .9];
        ylims.e3 = [6 20];
        ylims.e5 = [.5 .9];
    case 'rt-g'
        ylims.e0 = [.9 1.6];
        ylims.e3 = [.02 .08];
        ylims.e5 = [0 .7];
end

idx = [1 3 2];
nCV = numel(idx);

for iExp = 1:nExp
    expName = expNames{iExp};
    
    figure
    hold on
    b1 = bar(1:nCV, groupMean.(expName)(idx),'FaceColor',[.5 .5 .5]);
    p1 = errorbar(1:nCV, groupMean.(expName)(idx), groupSte.(expName)(idx),'k','LineStyle','none');
    ylim(ylims.(expName))
    title(expName)
end
