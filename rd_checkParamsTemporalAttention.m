% rd_checkParamsTemporalAttention.m

%% setup
% subject = 'ecPilot_cb_tilt1pt5_tc64-100_soa1000-1250';
subject = 'md_cbD10_tilt*_tc16-64_soa1000-1300';
runs = 1:3;

expName = 'E4_contrast_cbD10'; % 'E2_SOA_cbD6', 'E0_cb'
% dataDir = 'data';
% figDir = 'figures';
dataDir = pathToExpt('data');
figDir = pathToExpt('figures');
dataDir = sprintf('%s/%s/%s', dataDir, expName, subject(1:2));
figDir = sprintf('%s/%s/%s', figDir, expName, subject(1:2));

for iRun = 1:numel(runs)
    run = runs(iRun);
    %% load data
    dataFile = dir(sprintf('%s/%s_run%02d*TemporalAttention*', dataDir, subject, run));
    data = load(sprintf('%s/%s', dataDir, dataFile.name));
    
    p = data.expt.p;
    
    %% print params
    fprintf('\nRun %d', run)
    fprintf('\ntarget contrasts: %d %d\n', p.targetContrasts*100)
    fprintf('SOAs: %d %d\n', p.soas*1000)
    fprintf('target orientation: %1.1f %1.1f\n', p.targetOrientation)
    fprintf('cue validity factor: %d %d %d %d %d\n', p.cueValidityFactor)
end