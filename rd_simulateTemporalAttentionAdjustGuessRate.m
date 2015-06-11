% rd_simulateTemporalAttentionAdjustGuessRate.m

nSims = 10; % number of simulations

ns = [64 192]; % number of trials
sd = 10;
guessRates = [0.01 0.02 0.03 0.04 0.05 0.06];

model = Orientation(StandardMixtureModel, 2); % sd

% run simulation
for iN = 1:numel(ns)
    n = ns(iN);
    fprintf('\n\n[%s] n = %d\n', datestr(now), n)
    for iG = 1:numel(guessRates)
        g = guessRates(iG);
        fprintf('\n[%s] guess rate = %.02f\n', datestr(now), g)
        
        nG = round(n*g); % randomly guessing (contributing to the guess distribution)
        nSD = n-nG; % noisy estimation (contributing to the sd distribution)
        
        for iSim = 1:nSims
            fprintf('simulation %d\n', iSim)
            
            errors = [circ_vmrnd(0, deg2k(sd), nSD).*(180/pi); rand(nG,1)*180-90];
            
            fit = MemFit(errors, model, 'Verbosity', 0);
            
            gHat(iSim,iG) = fit.maxPosterior(1); % posteriorMean
            sdHat(iSim,iG) = fit.maxPosterior(2);
        end
    end
    results(iN).gHat = gHat;
    results(iN).sdHat = sdHat;
end

percentChange = (mean(results(1).gHat) - mean(results(2).gHat))./mean(results(2).gHat)*100;

figure
hold on
errorbar(guessRates, mean(results(1).gHat), std(results(1).gHat))
errorbar(guessRates, mean(results(2).gHat), std(results(2).gHat),'r')
plot([0 0.1], [0 0.1], 'k')
axis square
xlim([0 .1])
ylim([0 .1])
legend(num2str(ns'))
xlabel('true guess rate')
ylabel('estimated guess rate')
title('guess rate recovery for different numbers of trials')

% x = -90:90;
% y = (1-g).*vonmisespdf(x,mu,deg2k(sd))+(g).*1/180;
