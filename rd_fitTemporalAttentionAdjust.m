function [fit, err] = rd_fitTemporalAttentionAdjust(subjectID, run, saveData, bootRun)

% basic fitting
% data.errors = ?
% fit = MemFit(errors, model);
% fit = MemFit(data, model);

% orientation data
% MemFit(data, Orientation(WithBias(StandardMixtureModel), [1,3]))

% model comparison
% MemFit(data, {model1, model2})

if nargin < 4
    bootRun = 0;
end
if nargin < 3
    saveData = 0;
end

if bootRun > 0
    resample = 1;
else
    resample = 0;
end

separateConditions = 1; % 1 for normal, 0 for all conditions lumped together

%% setup
% subjectID = 'rd';
subject = sprintf('%s_a1_tc100_soa1000-1250', subjectID);
% run = 9;

expName = 'E3_adjust';
% dataDir = 'data';
% figDir = 'figures';
% dataDir = '~/Desktop/E3_data_rd_TOGO_20150821';
dataDir = pathToExpt('data');
figDir = pathToExpt('figures');
dataDir = sprintf('%s/%s/%s', dataDir, expName, subject(1:2));
figDir = sprintf('%s/%s/%s', figDir, expName, subject(1:2));

%% load data
dataFile = dir(sprintf('%s/%s_run%02d*', dataDir, subject, run));
load(sprintf('%s/%s', dataDir, dataFile(1).name))

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
if separateConditions
    for iEL = 1:2
        fprintf('\n%s', targetNames{iEL})
        for iV = 1:3
            fprintf('\n%s', validityNames{iV})
            errors = results.totals.all{iV,iEL}(:,errorIdx);
            
            if resample
                n = numel(errors);
                bootsam = ceil(n*rand(n,1));
                errors = errors(bootsam);
            end
            
            data.errors = errors';
            
            if ~isempty(strfind(modelName, 'swap'))
                [e, to, nto] = ...
                    rd_plotTemporalAttentionAdjustErrors(subjectID, run, 0);
                data.distractors = nto{iV,iEL}';
            end
            
            fit(iV,iEL) = MemFit(data, model, 'Verbosity', 0);
            %         fit(iV,iEL).mle = MLE(data, model);
            
            err{iV,iEL} = errors;
            
            %         PlotModelFit(model, fit(iV,iEL).mle, data, 'NewFigure', true);
        end
    end
else % lump data from all conditions together
    data.errors = [];
    if ~isempty(strfind(modelName, 'swap'))
        data.distractors = [];
    end
    for iEL = 1:2
        fprintf('\n%s', targetNames{iEL})
        for iV = 1:3
            fprintf('\n%s', validityNames{iV})
            errors = results.totals.all{iV,iEL}(:,errorIdx);
            
            if resample
                n = numel(errors);
                bootsam = ceil(n*rand(n,1));
                errors = errors(bootsam);
            end
            
            data.errors = [data.errors errors'];
            
            if ~isempty(strfind(modelName, 'swap'))
                [e, to, nto] = ...
                    rd_plotTemporalAttentionAdjustErrors(subjectID, run, 0);
                distractors = nto{iV,iEL}';
                data.distractors = [data.distractors distractors];
            end
        end
    end
    % fit all data at once
    fit = MemFit(data, model, 'Verbosity', 0);
    err = data.errors;
end

%% save data
if resample
    bootExt = sprintf('_boot%04d', bootRun);
    saveDir = sprintf('%s/bootstrap/%s', dataDir, modelName);
else
    bootExt = '';
    saveDir = dataDir;
end

if separateConditions==0
    lumpExt = '_lumped';
else
    lumpExt = '';
end

if saveData
    if ~exist(saveDir,'dir')
        mkdir(saveDir)
    end
    fileName = sprintf('%s_run%02d_%s%s%s.mat', subject, run, modelName, lumpExt, bootExt);
    save(sprintf('%s/%s',saveDir,fileName),'fit','err','model')
end


