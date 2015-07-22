% rd_simulateTemporalAttentionAdjustGuessRate2.m

nSims = 10; % number of simulations

ns = [64 192]; % number of trials
guessRates = [0.01 0.02 0.03 0.04 0.05 0.06];

%% specify model
modelName = 'variablePrecision';
switch modelName
    case 'fixedNoBias'
        model = Orientation(NoGuessingModel, 1); % sd
    case 'mixtureWithBias'
        model = Orientation(WithBias(StandardMixtureModel), [1,3]); % mu, sd
    case 'mixtureNoBias'
        model = Orientation(StandardMixtureModel, 2); % sd
        params = [0.05 10];
    case 'swapNoBias'
        model = Orientation(SwapModel,3); % sd
    case 'swapWithBias'
        model = Orientation(WithBias(SwapModel), [1 4]); % mu, sd
    case 'variablePrecision'
        model = Orientation(VariablePrecisionModel, [2,3]); % mnSTD, stdSTD
        params = [.05 10 3];
    case 'variablePrecisionGammaSD'
        model = Orientation(VariablePrecisionModel('HigherOrderDist','GammaSD'), [2,3]); % modeSTD, sdSTD
    case 'variablePrecisionNoGuess'
        model = Orientation(VariablePrecisionModel('HigherOrderDist','GaussianSDNoGuess'), [1,2]); % mnSTD, stdSTD    
    otherwise
        error('modelName not recognized')
end

%% run simulation
for iN = 1:numel(ns)
    n = ns(iN);
    fprintf('\n\n[%s] n = %d\n', datestr(now), n)
    for iG = 1:numel(guessRates)
        g = guessRates(iG);
        fprintf('\n[%s] guess rate = %.02f\n', datestr(now), g)
        params(1) = g;
        
        for iSim = 1:nSims
            fprintf('simulation %d\n', iSim)
            
            errors = SampleFromModel(model, params, [1 n]);
            
            fit(iSim,iG,iN) = MemFit(errors, model, 'Verbosity', 0);
        end
    end
end

%% organize results
paramNames = model.paramNames;
for iN = 1:numel(ns)
    for iG = 1:numel(guessRates)
        for iSim = 1:nSims
            for iP = 1:numel(params)
                paramName = paramNames{iP};
                results.(paramName)(iSim,iG,iN) = fit(iSim,iG,iN).maxPosterior(iP);
            end
        end
    end
end
                

%% plot
figure
hold on
errorbar(guessRates, mean(results.g(:,:,1)), std(results.g(:,:,1)))
errorbar(guessRates, mean(results.g(:,:,2)), std(results.g(:,:,2)), 'r')
plot([0 0.1], [0 0.1], 'k')
axis square
xlim([0 .1])
% ylim([0 .1])
legend(num2str(ns'))
xlabel('true guess rate')
ylabel('estimated guess rate')
title('guess rate recovery for different numbers of trials')

figure
hold on
errorbar(guessRates, mean(results.mnSTD(:,:,1)), std(results.mnSTD(:,:,1)))
errorbar(guessRates, mean(results.mnSTD(:,:,2)), std(results.mnSTD(:,:,2)), 'r')
plot([0 0.07], [10 10], 'k')
xlim([0 .07])
ylim([5 15])
legend(num2str(ns'))
xlabel('true guess rate')
ylabel('estimated mnSTD')

figure
hold on
errorbar(guessRates, mean(results.stdSTD(:,:,1)), std(results.stdSTD(:,:,1)))
errorbar(guessRates, mean(results.stdSTD(:,:,2)), std(results.stdSTD(:,:,2)), 'r')
plot([0 0.07], [3 3], 'k')
xlim([0 .07])
ylim([0 6])
legend(num2str(ns'))
xlabel('true guess rate')
ylabel('estimated stdSTD')
