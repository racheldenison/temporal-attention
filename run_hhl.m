% run_hhl.m

%% subject-specific settings
tilt = 2;
soas = [450 350 100 250 800 300 500 200 150 400];

%% what session?
fprintf('\nHello, Hsin-Hung!\n\n')

day = input('day? (1-10) ');
run = input('run? (0,1,2,3) ');
soa = soas(day);

%% message
fprintf('\nDay %d, run %d\nYour tilt is %d degrees, and your SOA is %d ms.\n\n', day, run, tilt, soa)
if run==0
    fprintf('This is your practice run, so stop after the\nfirst block (ctrl-C & clear all).\n\n')
end

ok = input('ok? (n if not) ','s');
if strcmp(ok,'n')
    error('Ok, try again!')
end

%% set parameters
p0.soas = [1000 1000+soa]/1000;
p0.targetOrientation = [-tilt tilt];
p0.targetStates = [-tilt tilt];

%% set runStr
if run==0
    runStr = 'PRAC2';
else
    runStr = sprintf('run%02d',run+3);
end

%% set subjectID
subjectID = sprintf('hl_cbD6_tilt%d_tc64_soa%d-%d_%s', tilt, round(p0.soas(1)*1000), round(p0.soas(2)*1000), runStr);

fprintf('Filename: %s\n\n', subjectID)

%% run experiment
% rd_temporalAttention(subjectID, [], p0)

