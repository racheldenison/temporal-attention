% rd_simulateTemporalAttentionAdjustFit.m

mu = 0;
sd = 25;
trials = [8 16 32 64 96 128 192 256];

model = Orientation(WithBias(StandardMixtureModel), [1,3]);

% run simulation
for iT = 1:numel(trials)
    nTrials = trials(iT);
    fprintf('[%s] %d trials\n', datestr(now), nTrials)
    
    errors = mu + sd.*randn(nTrials,1);
    
    fitParams(iT) = MemFit(errors, model, 'Verbosity', 2);
    fitGoodness(iT) = MemFit(errors, {model, model}, 'Verbosity', 2);
end
     
% organize results
for iT = 1:numel(trials)
    estParams(iT,:) = fitParams(iT).posteriorMean;
    estParamsLower(iT,:) = fitParams(iT).lowerCredible;
    estParamsUpper(iT,:) = fitParams(iT).upperCredible;
    estAIC(iT) = fitGoodness(iT).AIC(1);
    estLogLike(iT) = fitGoodness(iT).logLike(1);
end

% plots
for iP = 1:numel(model.paramNames)
    paramName = model.paramNames{iP};
    figure
    hold on
    if strcmp(paramName,'mu')
        plot([1 numel(trials)],[mu mu],'--r')
    elseif strcmp(paramName,'sd')
        plot([1 numel(trials)],[sd sd],'--r')
    end
    plot(estParams(:,iP),'k')
    plot(estParamsLower(:,iP),'Color',[.5 .5 .5])
    plot(estParamsUpper(:,iP),'Color',[.5 .5 .5])
    title(paramName)
    set(gca,'XTickLabel',trials)
    xlabel('number of trials')
end

figure
plot(estLogLike,'.-')
set(gca,'XTickLabel',trials)
xlabel('number of trials')
ylabel('log likelihood')

figure
plot(estAIC,'.-')
set(gca,'XTickLabel',trials)
xlabel('number of trials')
ylabel('AIC')