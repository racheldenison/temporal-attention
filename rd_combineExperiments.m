% rd_combineExperiments.m


%% setup
e0 = load('data/E0_workspace_run09_N10_20160224.mat');
e3 = load('data/E3_workspace_run09_N12_20160224.mat');
e5 = load('data/E5_workspace_run01_N10_20160806.mat');

for iT = 1:2
    data.e0(:,iT,:) = e0.accDataC{iT};
    data.e3 = e3.paramsData.sd*(-1);
    data.e5(:,iT,:) = e5.accDataC{iT};
end

expNames = {'e0','e3','e5'};
nExp = numel(expNames);

%% zscore across all data points per experiment
for iE = 1:nExp
    exp = expNames{iE};
    vals = data.(exp);
    m = mean(vals(:));
    sd = std(vals(:));
    zdata.(exp) = (vals-m)/sd;
end

%% make table
ztable = [];
for iE = 1:nExp
    exp = expNames{iE};
    vals = zdata.(exp);
    sz = size(vals);
    v = reshape(vals, prod(sz), 1);
    
    validity = repmat([1 3 2]', sz(2)*sz(3), 1);
    target = repmat([1 1 1 2 2 2]', sz(3), 1);
    subject = [1 1 1 1 1 1]'*(1:sz(3));
    experiment = ones(size(v))*iE;
    
    tab = [experiment subject(:) validity target v];
    ztable = [ztable; tab];
end

%% make table
dtable = [];
for iE = 1:nExp
    exp = expNames{iE};
    vals = data.(exp);
    sz = size(vals);
    v = reshape(vals, prod(sz), 1);
    
    validity = repmat([1 3 2]', sz(2)*sz(3), 1);
    target = repmat([1 1 1 2 2 2]', sz(3), 1);
    subject = [1 1 1 1 1 1]'*(1:sz(3));
    experiment = ones(size(v))*iE;
    
    tab = [experiment subject(:) validity target v];
    dtable = [dtable; tab];
end
    