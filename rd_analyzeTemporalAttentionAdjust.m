function [expt results] = rd_analyzeTemporalAttentionAdjust(expt, saveData, saveFigs, plotTimingFigs, saveTimingFigs, T1T2Axis, cleanRT)

if nargin < 7 || isempty(cleanRT)
    cleanRT = 0;
end
if nargin < 6 || isempty(T1T2Axis)
    T1T2Axis = 'all';
end
if nargin < 5 || isempty(saveTimingFigs)
    saveTimingFigs = 0;
end
if nargin < 4 || isempty(plotTimingFigs)
    plotTimingFigs = 0;
end
if nargin < 3 || isempty(saveFigs)
    saveFigs = 0;
end
if nargin < 2 || isempty(saveData)
    saveData = 0;
end

steOption = 'trial'; % 'trial','set'
fprintf('\nStandard error by %s\n\n', steOption)

%% Read out variables from expt
subjectID = expt.subjectID;
p = expt.p;
timing = expt.timing;
trials_headers = expt.trials_headers;
trials = expt.trials;
targetRotations = expt.targetRotations;

%% Get column idxs from trials_headers
targetContrastIdx = strcmp(trials_headers,'targetContrast');
respIntervalIdx = strcmp(trials_headers,'respInterval');
cueValidityIdx = strcmp(trials_headers,'cueValidity');
rtIdx = strcmp(trials_headers,'rt');
responseErrorIdx = strcmp(trials_headers,'responseError');
correctIdx = strcmp(trials_headers,'correct');

%% Clean RT if requested
if cleanRT
    rt0 = trials(:,rtIdx);
    cutoff = prctile(rt0,95);
    rt = rt0;
    rt(rt0 > cutoff) = NaN;
    trials(:,rtIdx) = rt;
    subjectID = [subjectID '_RTx'];
    
    % update expt
    expt.trials = trials;
    expt.subjectID = subjectID;
    
    figure
    subplot(2,1,1)
    hold on
    hist(rt0, 40)
    plot_vertical_line(cutoff);
    title('pre-clean')
    subplot(2,1,2)
    hold on
    hist(rt, 40)
    plot_vertical_line(cutoff);
    title('post-clean')
    xlabel('RT')
end

%% Description of trial sets
% *2 so we have 8 trials per rep for invalid/neutral
unitSet = ones(numel(p.targetContrasts)*numel(p.respInterval)*numel(p.cueValidityFactor)*2,1);
if strcmp(p.rotateTarget,'cb')
    setNum = unitSet*(1:size(trials,1)/4/size(unitSet,1));
    setNum = repmat(setNum(:),4,1);
else
    setNum = unitSet*(1:size(trials,1)/size(unitSet,1));
    setNum = setNum(:);
end
setNums = unique(setNum);

%% Selection of T1&T2 same/different/all axes
switch T1T2Axis
    case 'all'
        wAx = ones(size(targetRotations,1),1);
        axTitle = '';
    case 'same'
        wAx = abs(diff(targetRotations,1,2))<10;
        axTitle = 'T1 & T2 same axis';
    case 'diff'
        wAx = abs(diff(targetRotations,1,2))>10;
        axTitle = 'T1 & T2 different axes';
    otherwise
        error('T1T2Axis option not recognized')
end

%% Extra selection step if desired
% % only horizontal targets
% for i=1:size(targetRotations,1)
%     w0(i,1) = targetRotations(i,trials(i,2))==90;
% end

%% Analyze data
switch steOption
    case 'trial'
        %% standard error by trial
        for iRI = 1:numel(p.respInterval)
            for iCV = 1:numel(p.cueValidity)
                for iTC = 1:numel(p.targetContrasts)
                    w = wAx & trials(:,respIntervalIdx)==iRI & trials(:,cueValidityIdx)==iCV & trials(:,targetContrastIdx)==iTC;
                    
                    totals.all{iCV,iRI}(:,:,iTC) = trials(w,:);
                    
                    totals.means{iRI}(iCV,:,iTC) = nanmean(totals.all{iCV,iRI}(:,:,iTC),1);
                    totals.stds{iRI}(iCV,:,iTC) = nanstd(totals.all{iCV,iRI}(:,:,iTC),0,1);
                    totals.stes{iRI}(iCV,:,iTC) = totals.stds{iRI}(iCV,:,iTC)./sqrt(size(totals.all{iCV,iRI}(:,:,iTC),1));
                end
            end
        end
        
    case 'set'
        %% standard error by trial set
        for iSet = 1:numel(setNums)
            for iRI = 1:numel(p.respInterval)
                for iCV = 1:numel(p.cueValidity)
                    for iTC = 1:numel(p.targetContrasts)
                        w = setNum==iSet & wAx & trials(:,respIntervalIdx)==iRI & trials(:,cueValidityIdx)==iCV & trials(:,targetContrastIdx)==iTC;
                        
                        totals.all{iCV,iRI}(:,:,iTC,iSet) = trials(w,:);
                        
                        totals.setMeans{iRI}(iCV,:,iTC,iSet) = mean(totals.all{iCV,iRI}(:,:,iTC,iSet),1);
                    end
                end
                
                totals.means{iRI} = mean(totals.setMeans{iRI},4);
                totals.stds{iRI} = std(totals.setMeans{iRI},0,4);
                totals.stes{iRI} = totals.stds{iRI}./sqrt(numel(setNums));
            end
        end
    otherwise
        error('steOption not recognized')
end

%% Acc and RT means
for iRI = 1:numel(p.respInterval)
    errorMean{iRI} = squeeze(totals.means{iRI}(:,responseErrorIdx,:)); % [validity x contrast]
    errorSte{iRI} = squeeze(totals.stes{iRI}(:,responseErrorIdx,:));
    
    accMean{iRI} = squeeze(totals.means{iRI}(:,correctIdx,:)); % [validity x contrast]
    accSte{iRI} = squeeze(totals.stes{iRI}(:,correctIdx,:));
    
    rtMean{iRI} = squeeze(totals.means{iRI}(:,rtIdx,:));
    rtSte{iRI} = squeeze(totals.stes{iRI}(:,rtIdx,:));
end

%% Store data
results.totals = totals;
results.errorMean = errorMean;
results.errorSte = errorSte;
results.accMean = accMean;
results.accSte = accSte;
results.rtMean = rtMean;
results.rtSte = rtSte;
results.whenSaved = datestr(now);

%% Save data
if saveData
    fileName = sprintf('data/%s_TemporalAttentionAdjust_T1T2%s_%s.mat', subjectID, T1T2Axis, datestr(now, 'yyyymmdd'));
    save(fileName, 'expt', 'results')
end

%% Plot figs
intervalNames = {'early','late'};
accLims = [0 1];
rtLims = [0 10];
contrastLims = [p.targetContrasts(1)-0.05 p.targetContrasts(end)+0.05];
colors = get(0,'DefaultAxesColorOrder');

fig(1) = figure;
for iRI = 1:numel(p.respInterval)
    subplot(1,numel(p.respInterval),iRI)
    hold on
    plot(contrastLims, [0.5 0.5], '--k');
    
    if numel(p.targetContrasts)>1
        p1 = errorbar(repmat(p.targetContrasts',1,numel(p.cueValidity)),...
            errorMean{iRI}', errorSte{iRI}', '.', 'MarkerSize', 20);
    else
        for i = 1:length(errorMean{iRI})
            p1(i) = errorbar(p.targetContrasts,...
                errorMean{iRI}(i), errorSte{iRI}(i), '.', 'MarkerSize', 20);
            set(p1(i),'color', colors(i,:))
        end
    end
    xlabel('contrast')
    ylabel('error')
    legend(p1, num2str(p.cueValidity'),'location','best')
    title(intervalNames{iRI})
    xlim(contrastLims)
%     ylim(errorLims)
    rd_supertitle(subjectID);
    rd_raiseAxis(gca);
    rd_supertitle(axTitle);
end

fig(2) = figure;
for iRI = 1:numel(p.respInterval)
    subplot(1,numel(p.respInterval),iRI)
    hold on
    plot(contrastLims, [0.5 0.5], '--k');
    
    if numel(p.targetContrasts)>1
        p1 = errorbar(repmat(p.targetContrasts',1,numel(p.cueValidity)),...
            accMean{iRI}', accSte{iRI}', '.', 'MarkerSize', 20);
    else
        for i = 1:length(accMean{iRI})
            p1(i) = errorbar(p.targetContrasts,...
                accMean{iRI}(i), accSte{iRI}(i), '.', 'MarkerSize', 20);
            set(p1(i),'color', colors(i,:))
        end
    end
    xlabel('contrast')
    ylabel('acc')
    legend(p1, num2str(p.cueValidity'),'location','best')
    title(intervalNames{iRI})
    xlim(contrastLims)
    ylim(accLims)
    rd_supertitle(subjectID);
    rd_raiseAxis(gca);
    rd_supertitle(axTitle);
end

fig(3) = figure;
for iRI = 1:numel(p.respInterval)
    subplot(1,numel(p.respInterval),iRI)
    hold on
    
    if numel(p.targetContrasts)>1
        p1 = errorbar(repmat(p.targetContrasts',1,numel(p.cueValidity)),...
            rtMean{iRI}', rtSte{iRI}', '.', 'MarkerSize', 20);
    else
        for i = 1:length(rtMean{iRI})
            p1(i) = errorbar(p.targetContrasts,...
                rtMean{iRI}(i), rtSte{iRI}(i), '.', 'MarkerSize', 20);
            set(p1(i),'color', colors(i,:))
        end
    end

    xlabel('contrast')
    ylabel('rt')
    legend(num2str(p.cueValidity'),'location','best')
    title(intervalNames{iRI})
    xlim(contrastLims)
    ylim(rtLims)
    box off
    rd_supertitle(subjectID);
    rd_raiseAxis(gca);
    rd_supertitle(axTitle);
end

%% Save figs
if saveFigs
    figNames = {'error','acc','rt'};
    rd_saveAllFigs(fig, figNames, sprintf('%s_TemporalAttentionAdjust_T1T2%s', subjectID, T1T2Axis))
end

%% Plot timing
if plotTimingFigs
    tfig(1) = figure('Color','w');
    hold on
    nTrials = size(timing.dur.im1,1);
    plot(ones(nTrials,1), timing.dur.im1,'o')
    plot(2*ones(nTrials,1), timing.dur.im2,'o')
    plot(3*ones(nTrials,1), timing.dur.im1Im2SOA,'o')
    plot(4*ones(nTrials,1), timing.dur.cueIm1SOA,'o')
    plot(5*ones(nTrials,1), timing.dur.cueIm2SOA,'o')
    set(gca,'XTick',[1 2 3 4 5])
    set(gca,'XTickLabel',{'im 1','im 2','im1-im2 SOA','cue-im1 SOA','cue-im2 SOA'})
    ylabel('duration')
    
    tfig(2) = figure('Position',[1 1 700 300]);
    subplot(1,5,1)
    hist(timing.dur.im1)
    xlabel('im 1 duration (s)')
    ylabel('number of trials')
    subplot(1,5,2)
    hist(timing.dur.im2)
    xlabel('im 2 duration (s)')
    ylabel('number of trials')
    subplot(1,5,3)
    hist(timing.dur.im1Im2SOA)
    xlabel('im1-im2 SOA (s)')
    ylabel('number of trials')
    subplot(1,5,4)
    hist(timing.dur.cueIm1SOA)
    xlabel('cue-im1 SOA (s)')
    ylabel('number of trials')
    subplot(1,5,5)
    hist(timing.dur.cueIm2SOA)
    xlabel('cue-im2 SOA (s)')
    ylabel('number of trials')
    
%     tfig(3) = figure('Position',[1 1 700 300]);
%     fnames = fieldnames(timing.dur);
%     nFields = numel(fnames);
%     for iField = 1:nFields
%         subplot(nFields,1,iField)
%         plot(timing.dur.(fnames{iField}))
%         ylim([mean(timing.dur.(fnames{iField}))-0.01 mean(timing.dur.(fnames{iField}))+0.01])
%     end

    %% Save timing figs
    if saveTimingFigs
        figNames = {'timing','timingHist'};
        rd_saveAllFigs(tfig, figNames, sprintf('%s_TemporalAttentionAdjust', subjectID))
    end
end


