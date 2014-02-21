% rd_analyzeTemporalFields.m

subject = 'rdPilot1_t05_s20';
runs = 1:4;
nRuns = numel(runs);

for iRun = 1:nRuns
    run = runs(iRun);
    dataFile = dir(sprintf('data/%s*run%02d*', subject, run));
    load(sprintf('data/%s', dataFile.name))
    
    accData(:,iRun) = results.accMean;
end

soas = expt.p.soas;

accMean = mean(accData,2);
accSte = std(accData,0,2)./sqrt(nRuns);

figure
hold on
% plot(soas(1:end-1), accData(1:end-1,:))
errorbar(soas(1:end-1), accMean(1:end-1), accSte(1:end-1), '.-k')
errorbar(soas(end), accMean(end), accSte(end), '.-k')
xlabel('target-surround soa')
ylabel('accuracy')
title([subject ' runs ' num2str(runs)])

    
    