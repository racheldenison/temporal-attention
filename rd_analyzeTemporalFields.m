% rd_analyzeTemporalFields.m

subject = 'rdPilot1_t05_s20';
runs = 1:4;
nRuns = numel(runs);

for iRun = 1:nRuns
    run = runs(iRun);
    dataFile = dir(sprintf('data/%s_run%02d*', subject, run));
    load(sprintf('data/%s', dataFile.name))
    
    accData(:,iRun) = results.accMean;
    
    % results separated by target present/absent
    totals = results.totals;
    accPAData(:,1,iRun) = squeeze(mean(totals.all(1:2:end,8,:)));
    accPAData(:,2,iRun) = squeeze(mean(totals.all(2:2:end,8,:)));
 
end

soas = expt.p.soas;

accMean = mean(accData,2);
accSte = std(accData,0,2)./sqrt(nRuns);

accPAMean = mean(accPAData,3);
accPASte = std(accPAData,0,3)./sqrt(nRuns);

[dprime, criterion] = ccr_dprime(accPAMean(:,1), 1-accPAMean(:,2), 'yesno');

ylims = [0.3 1];

figure
hold on
% plot(soas(1:end-1), accData(1:end-1,:))
errorbar(soas(1:end-1), accMean(1:end-1), accSte(1:end-1), '.-k')
errorbar(soas(end), accMean(end), accSte(end), '.-k')
xlabel('target-surround soa')
ylabel('accuracy')
title([subject ' runs ' num2str(runs)])
ylim(ylims)

figure
hold on
errorbar(repmat(soas(1:end-1)',1,2), accPAMean(1:end-1,:), accPASte(1:end-1,:), '.-')
errorbar(soas(end), accPAMean(end,1), accPASte(end,1), '.-b')
errorbar(soas(end), accPAMean(end,2), accPASte(end,2), '.-g')
xlabel('target-surround soa')
ylabel('accuracy')
title([subject ' runs ' num2str(runs)])
legend('target present','target absent','location','best')
ylim(ylims)

figure
hold on
plot(soas(1:end-1), dprime(1:end-1), '.-k')
plot(soas(end), dprime(end), '.-k')
plot(soas(1:end-1), criterion(1:end-1), '.-b')
plot(soas(end), criterion(end), '.-b')
xlabel('target-surround soa')
ylabel('signal detection')
title([subject ' runs ' num2str(runs)])
legend('dprime','criterion')


    
    