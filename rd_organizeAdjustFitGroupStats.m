% rd_organizeAdjustFitGroupStats.m
%
% organize fit data for group stats

load data/adjust_fit_group_stats_run09_N10_20150216.mat

nSubjects = numel(subjectIDs);

targetNames = {'T1','T2'};
validityNames = {'valid','invalid','neutral'};

measures = fields(paramsData);

nVal = size(paramsData.sd,1);
nTarg = size(paramsData.sd,2);
nSub = size(paramsData.sd,3);

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

            