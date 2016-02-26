% rd_organizeAdjustFitGroupStats.m
%
% organize fit data for group stats

%% setup
% load data/adjust_fit_group_stats_run19_N10_20150216.mat
% load data/adjust_fit_group_stats_mixtureNoBias_run09_N10_20150303.mat
% load data/adjust_fit_group_stats_mixtureNoBias_run09_N12_20150407.mat
% load data/adjust_fit_group_stats_mixtureWithBiasMaxPosterior_run09_N12_20150512.mat
load data/adjust_fit_group_stats_mixtureNoBiasMaxPosterior_run09_N12_20150512.mat
% % load data/adjust_fit_group_stats_swapNoBiasMaxPosterior_run09_N12_20150512.mat
% load data/adjust_fit_group_stats_swapWithBiasMaxPosterior_run09_N12_20150516.mat
% load data/adust_fit_group_stats_variablePrecisionMaxPosterior_run09_N12_20150720.mat

% % load data/adjust_fit_group_stats_mixtureNoBiasMaxPosterior_run19_N12_20150513.mat
% load data/adjust_fit_group_stats_swapNoBiasMaxPosterior_run19_N12_20150601.mat
% load data/adjust_fit_group_stats_mixtureWithBiasMaxPosterior_run19_N12_20150601.mat

% load data/adjust_fit_group_stats_mixtureNoBiasMaxPosterior_run29_N12_20151209.mat

% paramsData paramsMean paramsSte run subjectIDs

subjects = [];

targetNames = {'T1','T2'};
validityNames = {'valid','invalid','neutral'};

% measures = fields(paramsData);
measures = {'g';'sd'};

if ~isempty(subjects)
    paramsData0 = paramsData;
    for iM = 1:numel(measures)
        m = measures{iM};
        paramsData.(m) = paramsData.(m)(:,:,subjects);
    end
end

nVal = size(paramsData.(measures{1}),1);
nTarg = size(paramsData.(measures{1}),2);
nSub = size(paramsData.(measures{1}),3);

%% make table
table_headers = [{'subject','T1T2','validity'}, measures'];
idx = 1;
for iSub = 1:nSub
    for iT = 1:nTarg
        for iVal = 1:nVal
            data = [];
            for iM = 1:numel(measures)
                m = measures{iM};
                data = [data paramsData.(m)(iVal,iT,iSub)];
            end
            table(idx,:) = [iSub iT iVal data];
            idx = idx + 1;
        end
    end
end

%% simple t-tests
for iM = 1:numel(measures)
    m = measures{iM};
    fprintf('\n\n%s', m)
    for iT = 1:nTarg
        dataV = squeeze(paramsData.(m)(1,iT,:));
        dataI = squeeze(paramsData.(m)(2,iT,:));
        dataN = squeeze(paramsData.(m)(3,iT,:));
    
        % qqplot to visualize deviations from normality
        figure
        subplot(1,3,1)
        qqplot(dataV-dataI)
        title('valid vs. invalid')
        subplot(1,3,2)
        qqplot(dataV-dataN)
        title('valid vs. neutral')
        subplot(1,3,3)
        qqplot(dataN-dataI)
        title('neutral vs. invalid')
        rd_supertitle(sprintf('%s %s', targetNames{iT}, m));
        rd_raiseAxis(gca);
        
        fprintf('\nT%d',iT)
        fprintf('\nt-test')
        [hVI pVI] = ttest(dataV,dataI);
        [hVN pVN] = ttest(dataV,dataN);
        [hNI pNI] = ttest(dataN,dataI);
        fprintf('\nvalid vs. invalid: p = %1.5f', pVI)
        fprintf('\nvalid vs. neutral: p = %1.5f', pVN)
        fprintf('\nneutral vs. invalid: p = %1.5f\n', pNI)
        
        fprintf('\nWilcoxon sign-rank test')
        [pVI] = signrank(dataV,dataI);
        [p, h, statsVI] = signrank(dataV,dataI,'method','approximate');
        [pVN] = signrank(dataV,dataN);
        [p, h, statsVN] = signrank(dataV,dataN,'method','approximate');
        [pNI] = signrank(dataN,dataI);
        [p, h, statsNI] = signrank(dataN,dataI,'method','approximate');
        fprintf('\nvalid vs. invalid: Z = %1.3f, p = %1.5f', statsVI.zval, pVI)
        fprintf('\nvalid vs. neutral: Z = %1.3f, p = %1.5f', statsVN.zval, pVN)
        fprintf('\nneutral vs. invalid: Z = %1.3f, p = %1.5f\n', statsNI.zval, pNI)
    end
end

%% Collapsing across T1 and T2
fprintf('\n\nCollapsing across T1 and T2\n')
for iM = 1:numel(measures)
    m = measures{iM};
    fprintf('\n%s\n', m)
    vals = squeeze(paramsData.(m)(:,1,:) + paramsData.(m)(:,2,:))./2;
    
    % qqplot to visualize deviations from normality
    figure
    subplot(1,3,1)
    qqplot(vals(1,:)-vals(2,:))
    title('valid vs. invalid')
    subplot(1,3,2)
    qqplot(vals(1,:)-vals(3,:))
    title('valid vs. neutral')
    subplot(1,3,3)
    qqplot(vals(2,:)-vals(3,:))
    title('neutral vs. invalid')
    rd_supertitle(sprintf('T1&T2 %s', m));
    rd_raiseAxis(gca);
    
    
    % distribution should be normal
    fprintf('t-test\n')
    [hvi pvi cvi svi] = ttest(vals(1,:),vals(2,:));
    [hvn pvn cvn svn] = ttest(vals(1,:),vals(3,:));
    [hni pni cni sni] = ttest(vals(2,:),vals(3,:));
    fprintf('valid vs. invalid, t(%d) = %1.3f, p = %1.4f\n', svi.df, svi.tstat, pvi)
    fprintf('valid vs. neutral, t(%d) = %1.3f, p = %1.4f\n', svn.df, svn.tstat, pvn)
    fprintf('neutral vs. invalid, t(%d) = %1.3f, p = %1.4f\n\n', sni.df, sni.tstat, pni)
    
    fprintf('Wilcoxon sign-rank test\n')
    % distribution should be symmetric about median
    [pvi hvi svi] = signrank(vals(1,:),vals(2,:)); % to return z-value, include 'method','approximate'
    [pvn hvn svn] = signrank(vals(1,:),vals(3,:));
    [pni hni sni] = signrank(vals(2,:),vals(3,:));
    fprintf('valid vs. invalid, p = %1.4f\n', pvi)
    fprintf('valid vs. neutral, p = %1.4f\n', pvn)
    fprintf('neutral vs. invalid, p = %1.4f\n\n', pni)
    
    fprintf('Sign test\n')
    [pvi] = signtest(vals(1,:),vals(2,:)); % to return z-value, include 'method','approximate'
    [pvn] = signtest(vals(1,:),vals(3,:));
    [pni] = signtest(vals(2,:),vals(3,:));
    fprintf('valid vs. invalid, p = %1.4f\n', pvi)
    fprintf('valid vs. neutral, p = %1.4f\n', pvn)
    fprintf('neutral vs. invalid, p = %1.4f\n\n', pni)
end

%% Ranomization tests
% load empirical null distribution
R = load('data/adjust_randomizationTest_workspace_run09_N12_20160108.mat');

% calculate observed pairwise differences
for iM = 1:numel(measures)
    m = measures{iM};
    pd(1,:,iM) = paramsMean.(m)(2,:) - paramsMean.(m)(1,:); % VI
    pd(2,:,iM) = paramsMean.(m)(3,:) - paramsMean.(m)(1,:); % VN
    pd(3,:,iM) = paramsMean.(m)(2,:) - paramsMean.(m)(3,:); % NI
end

fprintf('RANDOMIZATION TESTS: Fixed effects\n')
for iM = 1:numel(measures)
    m = measures{iM};
    fprintf('\n%s\n', m)
    for iT = 1:nTarg
        fprintf('T%d\n',iT)
        for iVC = 1:3 % validity comparison
            maxval = max(pd(iVC,iT,iM), -pd(iVC,iT,iM));
            minval = min(pd(iVC,iT,iM), -pd(iVC,iT,iM));
            pC(iVC,iT) = (nnz(R.pd(iVC,iT,iM,:) > maxval) + ...
                nnz(R.pd(iVC,iT,iM,:) < minval))/R.nSamples;
        end
        
        fprintf('valid vs. invalid, p = %1.3f\n', pC(1,iT))
        fprintf('valid vs. neutral, p = %1.3f\n', pC(2,iT))
        fprintf('neutral vs. invalid, p = %1.3f\n\n', pC(3,iT))
    end
end

%% randomization on subject means
nShuffles = 1000;
for iShuffle = 1:nShuffles
    for iSub = 1:nSub
        for iM = 1:numel(measures)
            m = measures{iM};
            newOrder = randperm(3);
            paramsDataShuffle.(m)(:,:,iSub,iShuffle) = ...
                paramsData.(m)(newOrder,:,iSub);
        end
    end
end

for iM = 1:numel(measures)
    m = measures{iM};
    paramsDataShuffleDiff.(m)(1,:,:,:) = paramsDataShuffle.(m)(2,:,:,:) - paramsDataShuffle.(m)(1,:,:,:); % VI
    paramsDataShuffleDiff.(m)(2,:,:,:) = paramsDataShuffle.(m)(3,:,:,:) - paramsDataShuffle.(m)(1,:,:,:); % VN
    paramsDataShuffleDiff.(m)(3,:,:,:) = paramsDataShuffle.(m)(2,:,:,:) - paramsDataShuffle.(m)(3,:,:,:); % NI
end

for iM = 1:numel(measures)
    m = measures{iM};
    for iVC = 1:3
        paramsMeanShuffleDiff.(m)(iVC,:,:) = squeeze(mean(paramsDataShuffleDiff.(m)(iVC,:,:,:),3));
    end
end

fprintf('RANDOMIZATION TESTS: Random effects\n')
for iM = 1:numel(measures)
    m = measures{iM};
    fprintf('\n%s\n', m)
    for iT = 1:nTarg
        fprintf('T%d\n',iT)
        for iVC = 1:3 % validity comparison
            maxval = max(pd(iVC,iT,iM), -pd(iVC,iT,iM));
            minval = min(pd(iVC,iT,iM), -pd(iVC,iT,iM));
            pC(iVC,iT) = (nnz(paramsMeanShuffleDiff.(m)(iVC,iT,:) > maxval) + ...
                nnz(paramsMeanShuffleDiff.(m)(iVC,iT,:) < minval))/nShuffles;
        end
        
        fprintf('valid vs. invalid, p = %1.3f\n', pC(1,iT))
        fprintf('valid vs. neutral, p = %1.3f\n', pC(2,iT))
        fprintf('neutral vs. invalid, p = %1.3f\n\n', pC(3,iT))
    end
end

%% Effect size
% calculate observed pairwise differences
m = 'sd';
for iT = 1:2
    pdData(1,:,iT) = paramsData.(m)(1,iT,:) - paramsData.(m)(2,iT,:); % VI
    pdData(2,:,iT) = paramsData.(m)(1,iT,:) - paramsData.(m)(3,iT,:); % VN
    pdData(3,:,iT) = paramsData.(m)(3,iT,:) - paramsData.(m)(2,iT,:); % NI
end
pdData = -pdData;

dP = mean(pdData,2)./std(pdData,0,2);

% R: pwr.t.test(d = 1.2335, sig.level = .05, power = .8, type = "paired")

