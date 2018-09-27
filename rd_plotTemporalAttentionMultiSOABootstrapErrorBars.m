% rd_plotTemporalAttentionMultiSOABootstrapErrorBars.m

%% Setup
subjectInit = 'jp';
exptName = 'cbD6'; % 'cbD6', 'cbD10'

run = 98;

% soa1 = [1000 1000 1000 1000 1000 1000];
% soa2 = [1150 1250 1400 1450 1500 1800];
soa1 = [1000 1000 1000 1000 1000 1000 1000 1000 1000];
soa2 = [1100 1200 1250 1300 1350 1400 1450 1500 1800];
soas = soa2 - soa1;
nSOA = numel(soas);

expName = 'E2_SOA_cbD6'; % 'E2_SOA_cbD6', 'E4_contrast_cbD10'
dataDir = pathToExpt('data');
dataDir = sprintf('%s/%s/%s', dataDir, expName, subjectInit(1:2));
figDir = sprintf('%s/%s/%s', pathToExpt('figures'), expName, subjectInit(1:2));

subjectID = sprintf('%s_%s*', subjectInit, exptName);

nboot = 1000;

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
end

%% Calculate confidence intervals
correctIdx = strcmp(expt.trials_headers,'correct');

for iSOA = 1:nSOA
    results = R{iSOA};
    for iT = 1:2
        for iV = 1:3
            vals = results.totals.all{iV,iT}(:,correctIdx);
            m = bootstrp(nboot, @mean, vals);
            ciBoot(iV,iT,iSOA,:) = prctile(m, [2.5 97.5]);
            accMean(iV,iT,iSOA) = results.accMean{iT}(iV);
            
            p = nnz(vals);
            n = numel(vals);
            b1 = betarnd(p,n-p+1,1000,1);
            b2 = betarnd(p+1,n-p,1000,1);
            ciBin(iV,iT,iSOA,1) = prctile(b1, 2.5);
            ciBin(iV,iT,iSOA,2) = prctile(b2, 97.5);
        end
    end
end

%% Plot
ci = ciBoot - repmat(accMean,[1,1,1,2]);
colors = {'b','r','k'};
figure
for iT = 1:2
    subplot(1,2,iT)
    hold on
    for iV = 1:3
        errorbar(soas, squeeze(accMean(iV,iT,:)), ...
            squeeze(ci(iV,iT,:,1)), squeeze(ci(iV,iT,:,2)), ...
            '.', 'Color', colors{iV});
    end
end





