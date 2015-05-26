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

% load data/adjust_fit_group_stats_mixtureNoBiasMaxPosterior_run19_N12_20150513.mat
% paramsData paramsMean paramsSte run subjectIDs

nSubjects = numel(subjectIDs);

targetNames = {'T1','T2'};
validityNames = {'valid','invalid','neutral'};

measures = fields(paramsData);

nVal = size(paramsData.sd,1);
nTarg = size(paramsData.sd,2);
nSub = size(paramsData.sd,3);

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
        [hVI pVI] = ttest(dataV,dataI);
        [hVN pVN] = ttest(dataV,dataN);
        [hNI pNI] = ttest(dataN,dataI);
        fprintf('\nT%d',iTarg)
        fprintf('\nvalid vs. invalid: p = %1.5f', pVI)
        fprintf('\nvalid vs. neutral: p = %1.5f', pVN)
        fprintf('\nneutral vs. invalid: p = %1.5f\n', pNI)
    end
end

            