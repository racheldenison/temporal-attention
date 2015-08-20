% rd_plotTemporalAttentionAdjustCompareFits.m

%% group i/o
subjectIDs = {'bl','rd','id','ec','ld','en','sj','ml','ca','jl','ew','jx'};
% subjectIDs = {'bl'};
run = 9;
nSubjects = numel(subjectIDs);

plotDistributions = 1;
saveFigs = 1;

groupFigTitle = [sprintf('%s ',subjectIDs{:}) sprintf('(N=%d), run %d', nSubjects, run)];

modelNames = {'mixtureNoBias','mixtureKurtosis'}; % 'fixedNoBias','mixtureWithBias','mixtureNoBias','swapNoBias', 'swapWithBias', 'variablePrecision', 'variablePrecisionGammaSD'

%% get data
for iSubject = 1:nSubjects
    %% indiv i/o
    subjectID = subjectIDs{iSubject};
    subject = sprintf('%s_a1_tc100_soa1000-1250', subjectID);

    expName = 'E3_adjust';
    % dataDir = 'data';
    % figDir = 'figures';
    dataDir = pathToExpt('data');
    figDir = pathToExpt('figures');
    dataDir = sprintf('%s/%s/%s', dataDir, expName, subject(1:2));
    figDir = sprintf('%s/%s/%s', figDir, expName, subject(1:2));
    
    %% load data
    dataFile = dir(sprintf('%s/%s_run%02d_%s_vs_%s.mat', dataDir, subject, run, modelNames{1}, modelNames{2}));
    load(sprintf('%s/%s', dataDir, dataFile.name))
    
    %% organize data
    measures = fields(fit);
    nMeasures = numel(measures);
    
    for iM = 1:nMeasures
        measure = measures{iM};
        groupData.(measure)(iSubject,:) = fit.(measure);
%         groupData.(measure) = real(groupData.(measure));
    end
end

%% summary stats
for iM = 1:nMeasures
    measure = measures{iM};
    groupMean(iM,:) = mean(groupData.(measure),1);
    groupSte(iM,:) = std(groupData.(measure),0,1)./sqrt(nSubjects);
    groupDiff(iM,:) = diff(groupData.(measure),1,2);
end

groupDiffMean = mean(groupDiff,2);
groupDiffSte = std(groupDiff,0,2)./sqrt(nSubjects);

%% plot
% individuals
for iM = 1:nMeasures
    measure = measures{iM};
    figure
    bar(groupData.(measure))
    title(measure)
    set(gca,'XTick',1:nSubjects)
    set(gca,'XTickLabel',subjectIDs)
    legend(modelNames)
end

% group
figure
hold on
barweb(groupMean, groupSte, [], measures)
legend(modelNames)

% individuals - differences between models
figure
for iM = 1:nMeasures
    measure = measures{iM};
    subplot(nMeasures,1,iM)
    bar(groupDiff(iM,:))
    title(measure)
    set(gca,'XTick',1:nSubjects)
    set(gca,'XTickLabel',subjectIDs)
end
rd_supertitle(sprintf('%s vs. %s (+ favors model 1)', modelNames{1}, modelNames{2}))
rd_raiseAxis(gca);


