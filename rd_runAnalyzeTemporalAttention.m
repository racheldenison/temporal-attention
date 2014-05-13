% rd_runAnalyzeTemporalAttention.m

subject = 'id_cb_tilt1pt5_tc64_soa1000-1250';
run = 9;

saveData = 0;
saveFigs = 1;
plotTimingFigs = 0;
saveTimingFigs = 0;

dataDir = pathToExpt('data');
% dataDir = 'data';
% dataDir = [pathToExpt('data') '/pilot/rd'];

% load data file
dataFile = dir(sprintf('%s/%s_run%02d*', ...
    dataDir, subject, run));

load(sprintf('%s/%s', dataDir, dataFile.name))

for t1t2 = {'same','diff'} % {'all'} 
	rd_analyzeTemporalAttention(expt, saveData, saveFigs, plotTimingFigs, saveTimingFigs, t1t2{1});
end
