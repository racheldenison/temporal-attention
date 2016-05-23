% rd_gatherTemporalAttentionMultiSOAData.m

%% Setup
subjectInit = 'ho';
exptName = 'cbD6'; % 'cbD6', 'cbD10'

run = 98;

soa1 = [1000 1000 1000 1000 1000 1000 1000 1000 1000 1000];
soa2 = [1100 1150 1200 1250 1300 1350 1400 1450 1500 1800];
soas = soa2 - soa1;
nSOA = numel(soas);

expName = 'E2_SOA_cbD6'; % 'E2_SOA_cbD6', 'E4_contrast_cbD10'
dataDir = pathToExpt('data');
dataDir = sprintf('%s/%s/%s', dataDir, expName, subjectInit(1:2));
figDir = sprintf('%s/%s/%s', pathToExpt('figures'), expName, subjectInit(1:2));

subjectID = sprintf('%s_%s*', subjectInit, exptName);

%% Get data
for iSOA = 1:numel(soa1)
    subject = sprintf('%s*_soa%d-%d', ...
        subjectID, soa1(iSOA), soa2(iSOA));
    
    % load data from a given soa
    dataFile = dir(sprintf('%s/%s_run%02d_T*.mat', ...
        dataDir, subject, run));
    if numel(dataFile)~=1
        fprintf('%s/%s_run%02d_T*', dataDir, subject, run)
        error('more or fewer than one matching data file')
    else
        load(sprintf('%s/%s', dataDir, dataFile.name))
    end
    
    R{iSOA} = results;
    E{iSOA} = expt;
end

%% Save
fileName = sprintf('data/%s_%s_run%d_alldata', expName, subjectInit, run);
save(fileName, 'R','E','soas','subjectInit','expName','run')

