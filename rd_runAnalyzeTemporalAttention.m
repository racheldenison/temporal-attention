% rd_runAnalyzeTemporalAttention.m

subject = 'idPilot_cb_tilt1pt5_tc64-100_soa1000-1250';
run = 2;

saveData = 0;
saveFigs = 1;
plotTimingFigs = 0;

% load data file
dataFile = dir(sprintf('%s/%s_run%02d*', ...
    pathToExpt('data'), subject, run));

load(sprintf('%s/%s',pathToExpt('data'), dataFile.name))

for t1t2 = {'all','same','diff'}
	rd_analyzeTemporalAttention(expt, saveData, saveFigs, plotTimingFigs, saveTimingFigs, t1t2{1});
end
