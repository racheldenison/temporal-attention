% rd_temporalAttentionMultiSOAStats.m

%% load data
D = load('data/E2_SOA_cbD6_run98_N4_workspace_20160128.mat');
% D = load('data/E2_SOA_cbD6_run98_N5_workspace_20180731.mat');
data = D.dpData;

subjects = D.subjectInits;
nSubjects = numel(subjects);

soas = D.t1t2soa;
nSOA = numel(soas);

ts = D.intervalNames;
nT = numel(ts);

vs = D.cueNames;
nV = numel(vs);

%% soa average
soaIdx = 3:6;
% soaIdx = 1:10;
dataSOAAve = [];
for iT = 1:nT
    dataSOAAve(:,:,iT) = mean(data{iT}(:,soaIdx,:),2);
end

%% pairwise comparisons
vals = dataSOAAve;
for iT = 1:nT    
    [hvi(iT), pvi(iT), civi(:,iT), statsvi(iT)] = ...
        ttest(vals(1,:,iT), vals(2,:,iT));
    
    [hvn(iT), pvn(iT), civn(:,iT), statsvn(iT)] = ...
        ttest(vals(1,:,iT), vals(3,:,iT));
    
    [hni(iT), pni(iT), cini(:,iT), statsni(iT)] = ...
        ttest(vals(3,:,iT), vals(2,:,iT));
end

fprintf('\nvalid vs. invalid: T1, p=%1.3f; T2, p=%1.3f', pvi)
fprintf('\nvalid vs. neutral: T1, p=%1.3f; T2, p=%1.3f', pvn)
fprintf('\nneutral vs. invalid: T1, p=%1.3f; T2, p=%1.3f\n\n', pni)

