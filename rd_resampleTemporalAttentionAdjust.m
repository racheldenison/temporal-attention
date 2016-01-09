function [fit, err] = rd_resampleTemporalAttentionAdjust(subjectID, run)

% [fit, err] = rd_resampleTemporalAttentionAdjust(subjectID, run)

%% setup
shuffleLabels = {'cueValidity'};

% subjectID = 'rd';
subject = sprintf('%s_a1_tc100_soa1000-1250', subjectID);
% run = 9;

expName = 'E3_adjust';
dataDir = pathToExpt('data');
figDir = pathToExpt('figures');
dataDir = sprintf('%s/%s/%s', dataDir, expName, subject(1:2));
figDir = sprintf('%s/%s/%s', figDir, expName, subject(1:2));

%% load data
dataFile = dir(sprintf('%s/%s_run%02d*', dataDir, subject, run));
load(sprintf('%s/%s', dataDir, dataFile(1).name))

%% shuffle labels and reanalyze
idx = zeros(size(shuffleLabels));
for iLabel = 1:numel(shuffleLabels)
    idx(iLabel) = find(strcmp(expt.trials_headers, shuffleLabels{iLabel}));
end
newOrder = randperm(size(expt.trials,1));
expt.trials(:,idx) = expt.trials(newOrder,idx);
[expt results] = rd_analyzeTemporalAttentionAdjust(expt);

%% specify model
modelName = 'mixtureNoBias';
switch modelName
    case 'fixedNoBias'
        model = Orientation(NoGuessingModel, 1); % sd
    case 'mixtureWithBias'
        model = Orientation(WithBias(StandardMixtureModel), [1,3]); % mu, sd
    case 'mixtureNoBias'
        model = Orientation(StandardMixtureModel, 2); % sd
    case 'swapNoBias'
        model = Orientation(SwapModel,3); % sd
    case 'swapWithBias'
        model = Orientation(WithBias(SwapModel), [1 4]); % mu, sd
    case 'variablePrecision'
        model = Orientation(VariablePrecisionModel, [2,3]); % mnSTD, stdSTD
    case 'variablePrecisionGammaSD'
        model = Orientation(VariablePrecisionModel('HigherOrderDist','GammaSD'), [2,3]); % modeSTD, sdSTD
    case 'variablePrecisionNoGuess'
        model = Orientation(VariablePrecisionModel('HigherOrderDist','GaussianSDNoGuess'), [1,2]); % mnSTD, stdSTD
    case 'mixtureKurtosis'
        model = Orientation(StandardMixtureModelVariableKurtosis, 2); % sd
    case 'mixtureKurtosisBySubject'
        allCondsFitFile = dir(sprintf('%s/%s_run%02d_mixtureKurtosis_lumped.mat', dataDir, subject, run));
        allConds = load(sprintf('%s/%s', dataDir, allCondsFitFile.name));
        n = allConds.fit.maxPosterior(strcmp(allConds.model.paramNames,'n'));
        fprintf('\n%s: n = %1.2f', modelName, n)
        model = Orientation(StandardMixtureModelFixedKurtosis(n), 2);
        model.n = n;
    otherwise
        error('modelName not recognized')
end

%% get errors and fit model
errorIdx = strcmp(expt.trials_headers, 'responseError');

targetNames = {'T1','T2'};
validityNames = {'valid','invalid','neutral'};
for iEL = 1:2
%     fprintf('\n%s', targetNames{iEL})
    for iV = 1:3
%         fprintf('\n%s', validityNames{iV})
        errors = results.totals.all{iV,iEL}(:,errorIdx);
        
        data.errors = errors';
        
        if ~isempty(strfind(modelName, 'swap'))
            [e, to, nto] = ...
                rd_plotTemporalAttentionAdjustErrors(subjectID, run, 0);
            data.distractors = nto{iV,iEL}';
        end
        
%         fit(iV,iEL) = MemFit(data, model, 'Verbosity', 0);
        fit(iV,iEL).mle = MLE(data, model);
        
        err{iV,iEL} = errors;
        
%         PlotModelFit(model, fit(iV,iEL).mle, data, 'NewFigure', true);
    end
end


