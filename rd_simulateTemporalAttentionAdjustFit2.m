% rd_simulateTemporalAttentionAdjustFit2.m

% trials = [64 96 128 192 256 320 384];
trials = [64 192];

%% specify model
modelName = 'mixtureNoBias';
switch modelName
    case 'fixedNoBias'
        model = Orientation(NoGuessingModel, 1); % sd
    case 'mixtureWithBias'
        model = Orientation(WithBias(StandardMixtureModel), [1,3]); % mu, sd
    case 'mixtureNoBias'
        model = Orientation(StandardMixtureModel, 2); % sd
        params = [0.05 15];
    case 'swapNoBias'
        model = Orientation(SwapModel,3); % sd
    case 'swapWithBias'
        model = Orientation(WithBias(SwapModel), [1 4]); % mu, sd
    case 'variablePrecision'
        model = Orientation(VariablePrecisionModel, [2,3]); % mnSTD, stdSTD
        params = [.05 15 3];
    case 'variablePrecisionGammaSD'
        model = Orientation(VariablePrecisionModel('HigherOrderDist','GammaSD'), [2,3]); % modeSTD, sdSTD
    case 'variablePrecisionNoGuess'
        model = Orientation(VariablePrecisionModel('HigherOrderDist','GaussianSDNoGuess'), [1,2]); % mnSTD, stdSTD    
    otherwise
        error('modelName not recognized')
end

%% run simulation
for iT = 1:numel(trials)
    nTrials = trials(iT);
    fprintf('[%s] %d trials\n', datestr(now), nTrials)
    errors = SampleFromModel(model, params, [1 nTrials]);
    % figure, hist(errors,-179:180)
    fit(iT) = MemFit(errors, model,'Verbosity',2);
end

%% organize results
for iT = 1:numel(trials)
    estParams(iT,:) = fit(iT).posteriorMean;
    estParamsLower(iT,:) = fit(iT).lowerCredible;
    estParamsUpper(iT,:) = fit(iT).upperCredible;
end

%% plots
for iP = 1:numel(model.paramNames)
    paramName = model.paramNames{iP};
    figure
    hold on
    plot([trials(1) trials(end)],[params(iP) params(iP)],'--r')
    plot(trials, estParams(:,iP),'k')
    plot(trials, estParamsLower(:,iP),'Color',[.5 .5 .5])
    plot(trials, estParamsUpper(:,iP),'Color',[.5 .5 .5])
    title(paramName)
    set(gca,'XTick',trials)
    set(gca,'XTickLabel',trials)
    xlabel('number of trials')
end

%% something else
% trials = 192; guessRate = .1;
% guessTrials = round(trials*guessRate); 
% ngTrials = trials - guessTrials; 
% vpData = [SampleFromModel(VariablePrecisionModel,{0,15,3},[1 ngTrials]) ...
%     randsample(-89:90,guessTrials,'true')];

%% plot figs
% figure
% subplot(1,2,1)
% hist(errors)
% subplot(1,2,2)
% hist(errors(abs(errors)<100))


