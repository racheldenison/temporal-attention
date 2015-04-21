function [fit, err] = rd_fitTemporalAttentionAdjust(subjectID, run, saveData)

% basic fitting
% data.errors = ?
% fit = MemFit(errors, model);
% fit = MemFit(data, model);

% orientation data
% MemFit(data, Orientation(WithBias(StandardMixtureModel), [1,3]))

% model comparison
% MemFit(data, {model1, model2})

if nargin < 3
    saveData = 0;
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
modelName = 'mixtureWithBias';
switch modelName
    case 'mixtureWithBias'
        model = Orientation(WithBias(StandardMixtureModel), [1,3]); % mu, sd
    case 'mixtureNoBias'
        model = Orientation(StandardMixtureModel, 2); % sd
    case 'variablePrecision'
        model = Orientation(VariablePrecisionModel, [2,3]); % mnSTD, stdSTD
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
        
        fit(iV,iEL) = MemFit(errors, model, 'Verbosity', 0);
        
        err{iV,iEL} = errors;
    end
end

%% save data
if saveData
    fileName = sprintf('%s_run%02d_%s.mat', subject, run, modelName);
    save(sprintf('%s/%s',dataDir,fileName),'fit','err','model')
end


