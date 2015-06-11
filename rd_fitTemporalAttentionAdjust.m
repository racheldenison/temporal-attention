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

%% setup
% subjectID = 'rd';
subject = sprintf('%s_a1_tc100_soa1000-1250', subjectID);
% run = 9;

expName = 'E3_adjust';
% dataDir = 'data';
% figDir = 'figures';
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
    case 'mixtureWithBias'
        model = Orientation(WithBias(StandardMixtureModel), [1,3]); % mu, sd
    case 'mixtureNoBias'
        model = Orientation(StandardMixtureModel, 2); % sd
    case 'variablePrecision'
        model = Orientation(VariablePrecisionModel, [2,3]); % mnSTD, stdSTD
    case 'swapNoBias'
        model = Orientation(SwapModel,3); % sd
    case 'swapWithBias'
        model = Orientation(WithBias(SwapModel), [1 4]); % mu, sd
    otherwise
        error('modelName not recognized')
end

%% get errors and fit model
errorIdx = strcmp(expt.trials_headers, 'responseError');

targetNames = {'T1','T2'};
validityNames = {'valid','invalid','neutral'};
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
        
        err{iV,iEL} = errors;
    end
end

%% save data
if resample
    bootExt = sprintf('_boot%04d', bootRun);
    saveDir = sprintf('%s/bootstrap/%s', dataDir, modelName);
else
    bootExt = '';
    saveDir = dataDir;
end

if saveData
    if ~exist(saveDir,'dir')
        mkdir(saveDir)
    end
    fileName = sprintf('%s_run%02d_%s%s.mat', subject, run, modelName, bootExt);
    save(sprintf('%s/%s',saveDir,fileName),'fit','err','model')
end


