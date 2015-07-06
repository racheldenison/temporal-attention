function [fit, err] = rd_fitTemporalAttentionAdjustVP(subjectID, run, saveData, bootRun)

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

plotFigs = 0;

warning('off', 'MATLAB:divideByZero')

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

errorIdx = strcmp(expt.trials_headers, 'responseError');

%% specify model
modelName = 'VPK';
switch modelName
    case 'VP'
        fixedKappa = 0;
    case 'VPK'
        fixedKappa = 1;
        kappa_r = 4.5139;
        if exist('kappa_r','var')
            data.kappa_r = kappa_r;
        else
            % combine all errors to fix kappa_r
            allErrors = [];
            for iEL = 1:2
                for iV = 1:3
                    allErrors = [allErrors; results.totals.all{iV,iEL}(:,errorIdx)];
                end
            end
            allData.error_vec = 2*(allErrors*pi/180)';
            allData.N = 2*ones(size(allData.error_vec));
            for iFit = 1:100
                fprintf('%d ', iFit)
                allfitpars(iFit,:) = fit_VPA_model_rd(allData);
            end
            data.kappa_r = median(allfitpars(:,4));
        end
    otherwise
        error('modelName not recognized')
end

%% get errors and fit model
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
        
        data.error_vec = 2*(errors*pi/180)';
        data.N = 2*ones(size(data.error_vec));
        
%         if ~isempty(strfind(modelName, 'swap'))
%             [e, to, nto] = ...
%                 rd_plotTemporalAttentionAdjustErrors(subjectID, run, 0);
%             data.distractors = nto{iV,iEL}';
%         end
        
        [fitpars, max_lh, AIC, BIC] = fit_VPA_model_rd(data);
        
        fit(iV,iEL).params = fitpars;
        fit(iV,iEL).maxLH = max_lh;
        fit(iV,iEL).AIC = AIC;
        fit(iV,iEL).BIC = BIC;
        
        err{iV,iEL} = data.error_vec;
    end
end

%% view sample fit
if plotFigs
    v = 3; el = 2;
    fitpars = fit(v,el).params;
    error_vec = err{v,el};
    
    data_fit = gen_fake_VPA_data(fitpars,1e5,2);
    
    % plot fit
    figure
    X = linspace(-pi,pi,52);
    X = X(1:end-1)+diff(X(1:2))/2;
    Y_emp = hist(error_vec,X);
    Y_emp = Y_emp/sum(Y_emp)/diff(X(1:2));
    Y_fit = hist(data_fit.error_vec,X);
    Y_fit = Y_fit/sum(Y_fit)/diff(X(1:2));
    bar(X,Y_emp,'k');
    hold on
    plot(X,Y_fit,'r-','Linewidth',3)
    legend('Data','Fit')
    xlabel('Response error');
    ylabel('Probability');
    xlim([-pi pi]);
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
    save(sprintf('%s/%s',saveDir,fileName),'fit','err')
end


