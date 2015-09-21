% rd_organizeAdjustFitGroupStats.m
%
% organize fit data for group stats

%% setup
% load data/adjust_fit_group_stats_run19_N10_20150216.mat
% load data/adjust_fit_group_stats_mixtureNoBias_run09_N10_20150303.mat
% load data/adjust_fit_group_stats_mixtureNoBias_run09_N12_20150407.mat
% load data/adjust_fit_group_stats_mixtureWithBiasMaxPosterior_run09_N12_20150512.mat
load data/adjust_fit_group_stats_mixtureNoBiasMaxPosterior_run09_N12_20150512.mat
% load data/adjust_fit_group_stats_swapNoBiasMaxPosterior_run09_N12_20150512.mat
% load data/adjust_fit_group_stats_swapWithBiasMaxPosterior_run09_N12_20150516.mat
% load data/adust_fit_group_stats_variablePrecisionMaxPosterior_run09_N12_20150720.mat

% load data/adjust_fit_group_stats_mixtureNoBiasMaxPosterior_run19_N12_20150513.mat
% load data/adjust_fit_group_stats_swapNoBiasMaxPosterior_run19_N12_20150601.mat
% load data/adjust_fit_group_stats_mixtureWithBiasMaxPosterior_run19_N12_20150601.mat
% paramsData paramsMean paramsSte run subjectIDs

nSubjects = numel(subjectIDs);

targetNames = {'T1','T2'};
validityNames = {'valid','invalid','neutral'};

measures = fields(paramsData);

nVal = size(paramsData.(measures{1}),1);
nTarg = size(paramsData.(measures{1}),2);
nSub = size(paramsData.(measures{1}),3);

%% make table
table_headers = [{'subject','T1T2','validity'}, measures'];
idx = 1;
for iSub = 1:nSub
    for iTarg = 1:nTarg
        for iVal = 1:nVal
            data = [];
            for iM = 1:numel(measures)
                m = measures{iM};
                data = [data paramsData.(m)(iVal,iTarg,iSub)];
            end
            table(idx,:) = [iSub iTarg iVal data];
            idx = idx + 1;
        end
    end
end

%% simple t-tests
for iM = 1:numel(measures)
    m = measures{iM};
    fprintf('\n\n%s', m)
    for iTarg = 1:nTarg
        dataV = squeeze(paramsData.(m)(1,iTarg,:));
        dataI = squeeze(paramsData.(m)(2,iTarg,:));
        dataN = squeeze(paramsData.(m)(3,iTarg,:));
    
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
        rd_supertitle(sprintf('%s %s', targetNames{iTarg}, m));
        rd_raiseAxis(gca);
        
        fprintf('\nT%d',iTarg)
        fprintf('\nt-test')
        [hVI pVI] = ttest(dataV,dataI);
        [hVN pVN] = ttest(dataV,dataN);
        [hNI pNI] = ttest(dataN,dataI);
        fprintf('\nvalid vs. invalid: p = %1.5f', pVI)
        fprintf('\nvalid vs. neutral: p = %1.5f', pVN)
        fprintf('\nneutral vs. invalid: p = %1.5f\n', pNI)
        
        fprintf('\nWilcoxon sign-rank test')
        [pVI] = signrank(dataV,dataI);
        [pVN] = signrank(dataV,dataN);
        [pNI] = signrank(dataN,dataI);
        fprintf('\nvalid vs. invalid: p = %1.5f', pVI)
        fprintf('\nvalid vs. neutral: p = %1.5f', pVN)
        fprintf('\nneutral vs. invalid: p = %1.5f\n', pNI)
    end
end

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



            