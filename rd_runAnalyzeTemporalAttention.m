% rd_runAnalyzeTemporalAttention.m

subject = 'rdPilot_cb_tilt2_tc10-20_soa1000-1250';
run = 2;

saveData = 0;
saveFigs = 1;

% load data file
dataFile = dir(sprintf('%s/%s_run%02d*', ...
    pathToExpt('data'), subject, run));

load(sprintf('%s/%s',pathToExpt('data'), dataFile.name))

for t1t2 = {'same','diff'}
	rd_analyzeTemporalAttention(expt, saveData, saveFigs, t1t2{1});
end
