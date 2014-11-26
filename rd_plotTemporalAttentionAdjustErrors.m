function [errors, targetOrients, nonTargetOrients, targetOrientDiff] = ...
    rd_analyzeTemporalAttentionAdjustErrors(subjectID)

%% setup
% subjectID = 'rd';
subject = sprintf('%s_a1_tc100_soa1000-1250', subjectID);
run = 9;

expName = 'E3_adjust';
% dataDir = 'data';
% figDir = 'figures';
dataDir = pathToExpt('data');
figDir = pathToExpt('figures');
dataDir = sprintf('%s/%s/%s', dataDir, expName, subject(1:2));
figDir = sprintf('%s/%s/%s', figDir, expName, subject(1:2));

%% load data
dataFile = dir(sprintf('%s/%s_run%02d*', dataDir, subject, run));
load(sprintf('%s/%s', dataDir, dataFile.name))

%% get trials info
p = expt.p;
trials = expt.trials;
trials_headers = expt.trials_headers;
targetRotations = expt.targetRotations;

targetContrastIdx = strcmp(trials_headers,'targetContrast');
respIntervalIdx = strcmp(trials_headers,'respInterval');
cueValidityIdx = strcmp(trials_headers,'cueValidity');
responseErrorIdx = strcmp(trials_headers,'responseError');

%% group trials and target rotations by validity
for iRI = 1:numel(p.respInterval)
    for iCV = 1:numel(p.cueValidity)
        for iTC = 1:numel(p.targetContrasts)
            w = trials(:,respIntervalIdx)==iRI & trials(:,cueValidityIdx)==iCV & trials(:,targetContrastIdx)==iTC;
            
            totals.all{iCV,iRI}(:,:,iTC) = trials(w,:);
            totals.targetRot{iCV,iRI}(:,:,iTC) = targetRotations(w,:);
        end
    end
end

%% analyze how errors depend on the non-post-cued target
for iRI = 1:numel(p.respInterval)
    for iCV = 1:numel(p.cueValidity)
        rotations = totals.targetRot{iCV,iRI};
        
        targetIdx = totals.all{iCV,iRI}(1,respIntervalIdx); % taking advantage of the previous grouping
        targetOrients{iCV,iRI} = rotations(:,targetIdx);
        nonTargetOrients{iCV,iRI} = rotations(:,3-targetIdx);
        errors{iCV,iRI} = totals.all{iCV,iRI}(:,responseErrorIdx);
        
        tod = nonTargetOrients{iCV,iRI} - targetOrients{iCV,iRI};
        
        tod(tod>90) = tod(tod>90) - 180;
        tod(tod<-90) = tod(tod<-90) + 180;
        
        targetOrientDiff{iCV,iRI} = tod;
    end
end

%% plot figs
targetNames = {'T1','T2'};
colors = {'b','g','r'};
errorLims = [-100 100];
orientationLims = [-10 190];
errorXTicks = [-90 -45 0 45 90];
orientationXTicks = [0 45 90 135 180];

% target orientation
figure
for iRI = 1:numel(p.respInterval)
    subplot(numel(p.respInterval),1,iRI)
    hold on
    plot(orientationLims,[0 0], 'k')
    
    for iCV = 1:numel(p.cueValidity)
        plot(targetOrients{iCV,iRI}, errors{iCV,iRI}, '.', 'Color', colors{iCV})
    end
    
    set(gca,'XTick', orientationXTicks)
    xlim(orientationLims)
    ylim(errorLims)
    ylabel('error')
    title(targetNames{iRI})
end
xlabel('target orientation')
rd_supertitle(sprintf('%s run %d', subjectID, run));
rd_raiseAxis(gca);

% orientation difference between targets
figure
for iRI = 1:numel(p.respInterval)
    subplot(numel(p.respInterval),1,iRI)
    hold on
    plot(errorLims,[0 0], 'k')
    
    for iCV = 1:numel(p.cueValidity)
        plot(targetOrientDiff{iCV,iRI}, errors{iCV,iRI}, '.', 'Color', colors{iCV})
    end
    
    set(gca,'XTick', errorXTicks)
    xlim(errorLims)
    ylim(errorLims)
    ylabel('error')
    title(targetNames{iRI})
end
xlabel('non-target - target orientation difference')
rd_supertitle(sprintf('%s run %d', subjectID, run));
rd_raiseAxis(gca);


