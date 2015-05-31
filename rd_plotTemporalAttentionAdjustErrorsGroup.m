% rd_plotTemporalAttentionAdjustErrorsGroup.m

%% setup
subjectIDs = {'bl','rd','id','ec','ld','en','sj','ml','ca','jl','ew','jx'};
nSubjects = numel(subjectIDs);

run = 9;
plotIndivFigs = 0;

analyzeProbe = 0;

%% get data
for iSubject = 1:nSubjects
    subjectID = subjectIDs{iSubject};
    [groupData0(iSubject).errors, ...
        groupData0(iSubject).targetOrients, ...
        groupData0(iSubject).nonTargetOrients, ...
        groupData0(iSubject).targetOrientDiff, ...
        groupData0(iSubject).probeOrients, ...
        groupData0(iSubject).probeOrientDiff, ...
        groupData0(iSubject).responses, ...
        groupData0(iSubject).targetOrientDiffSmooth] = ...
            rd_plotTemporalAttentionAdjustErrors(subjectID, run, plotIndivFigs);
end

if analyzeProbe==0
    groupData0 = rmfield(groupData0,'probeOrients');
    groupData0 = rmfield(groupData0,'probeOrientDiff');
end

%% organize data
fNames = fieldnames(groupData0(1));
for iSubject = 1:nSubjects
    for iF = 1:numel(fNames)
        for iRI = 1:2
            for iCV = 1:3
                fName = fNames{iF};
                switch fName
                    case 'targetOrientDiffSmooth'
                        groupData.(fName)(:,:,:,iSubject) = groupData0(iSubject).(fName);
                    otherwise
                        groupData.(fName){iCV,iRI}(:,iSubject) = ...
                            groupData0(iSubject).(fName){iCV,iRI};
                end
            end
        end
    end
end

% calculate mean errors for each x-axis value
for iRI = 1:2
    for iCV = 1:3
        allErrors = groupData.errors{iCV,iRI}(:);
        allTargetOrients = groupData.targetOrients{iCV,iRI}(:);
        allNonTargetOrients = groupData.nonTargetOrients{iCV,iRI}(:);
        allTargetOrientDiff = groupData.targetOrientDiff{iCV,iRI}(:);
        if analyzeProbe
            allProbeOrients = groupData.probeOrients{iCV,iRI}(:);
            allProbeOrientDiff = groupData.probeOrientDiff{iCV,iRI}(:);
        end
        
        ods = unique(allTargetOrientDiff);
        for iOD = 1:numel(ods)
            odIdx = allTargetOrientDiff==ods(iOD);
            errorByOD{iCV,iRI}(iOD) = mean(allErrors(odIdx));
            errorByODStd{iCV,iRI}(iOD) = std(allErrors(odIdx));
            errorByODAbs{iCV,iRI}(iOD) = mean(abs(allErrors(odIdx)));
        end
        
        tos = unique(allTargetOrients);
        for iTO = 1:numel(tos)
            toIdx = allTargetOrients==tos(iTO);
            errorByTO{iCV,iRI}(iTO) = mean(allErrors(toIdx));
        end
        
        ntos = unique(allNonTargetOrients);
        for iNTO = 1:numel(ntos)
            ntoIdx = allNonTargetOrients==ntos(iNTO);
            errorByNTO{iCV,iRI}(iNTO) = mean(allErrors(ntoIdx));
        end
        
        % joint of target and non-target
        binSize = 10;
        orients = 0:binSize:179;
        for iTO = 1:numel(orients)
            to = orients(iTO);
            for iNTO = 1:numel(orients)
                nto = orients(iNTO);
                idxL = allTargetOrients>=to & allNonTargetOrients>=nto;
                idxU = allTargetOrients<to+binSize & allNonTargetOrients<nto+binSize;
                idx = idxL & idxU;
                if nnz(idx)==0
                    val = NaN;
                else
                    val = mean(allErrors(idx));
                end
                errorByTNTO{iCV,iRI}(iTO,iNTO) = val;
            end
        end       
        
        if analyzeProbe
            pos = unique(allProbeOrients);
            for iPO = 1:numel(pos)
                poIdx = allProbeOrients==pos(iPO);
                errorByPO{iCV,iRI}(iPO) = mean(allErrors(poIdx));
            end
            
            pds = unique(allProbeOrientDiff);
            for iPD = 1:numel(pds)
                pdIdx = allProbeOrientDiff==pds(iPD);
                errorByPD{iCV,iRI}(iPD) = mean(allErrors(pdIdx));
            end
        else
            pos = NaN;
            pds = NaN;
        end
        
        odx{iCV,iRI} = ods;
        tox{iCV,iRI} = tos;
        ntox{iCV,iRI} = ntos;
        pox{iCV,iRI} = pos;
        pdx{iCV,iRI} = pds;
    end
end

% weighted average of target x non-target error matrices
allErrorByTNTO = (errorByTNTO{1,1} + errorByTNTO{1,2}).*0.6 + ...
    (errorByTNTO{2,1} + errorByTNTO{2,2}).*0.2 + ...
    (errorByTNTO{3,1} + errorByTNTO{3,2}).*0.2;

%% smoothed data
groupMean.targetOrientDiffSmooth = nanmean(groupData.targetOrientDiffSmooth, 4);
groupSte.targetOrientDiffSmooth = nanstd(groupData.targetOrientDiffSmooth, 0, 4)./sqrt(nSubjects);

winSize = 29; % check in rd_plotTemporalAttentionAdjustErrors.m
steps = -90+floor(winSize/2):90-floor(winSize/2);

%% plot figures
targetNames = {'T1','T2'};
colors = {'b','g','r'};
errorLims = [-100 100];
orientationLims = [-10 190];
errorXTicks = [-90 -45 0 45 90];
orientationXTicks = [0 45 90 135 180];

smoothSize = 10; % 5
b = (1/smoothSize)*ones(1,smoothSize);
a = 1;

% target orientation
figure
for iRI = 1:2
    subplot(2,1,iRI)
    hold on
    plot(orientationLims,[0 0], 'k')
    
    for iCV = 1:3
%         w = groupData.targetOrients{iCV,iRI}>90;
        w = logical(ones(size(groupData.targetOrients{iCV,iRI})));
        plot(groupData.targetOrients{iCV,iRI}(w), groupData.errors{iCV,iRI}(w), '.', 'Color', colors{iCV})
    end
    
    for iCV = 1:3
        smoothError = filter(b,a,errorByTO{iCV,iRI});
        plot(tox{iCV,iRI}, smoothError,'Color', colors{iCV}, 'LineWidth',2)
    end
    
    set(gca,'XTick', orientationXTicks)
    xlim(orientationLims)
    ylim(errorLims)
    ylabel('error')
    title(targetNames{iRI})
end
xlabel('target orientation')
rd_supertitle(sprintf('%s ', subjectIDs{:}));
rd_raiseAxis(gca);

% flip flip
figure
for iRI = 1:2
    subplot(2,1,iRI)
    hold on
    for iCV = 1:3
        to = groupData.targetOrients{iCV,iRI};
        e = groupData.errors{iCV,iRI};
        w = to>90;
        targetOrientsFF{iCV,iRI} = to;
        targetOrientsFF{iCV,iRI}(w) = 180 - to(w);
        errorsFF{iCV,iRI} = e;
        errorsFF{iCV,iRI}(w) = -e(w);
        
        eff = errorsFF{iCV,iRI}(:);
        toff = targetOrientsFF{iCV,iRI}(:);
        
        plot(toff, eff, '.', 'Color', colors{iCV})
    end
    plot([0 90], [0 0], 'k', 'LineWidth', 2)
    xlim([0 90])
    xlabel('flipflip target orientation')
    ylabel('flipflip error')
end
toErrorsFF = errorsFF;


% non-target orientation
figure
for iRI = 1:2
    subplot(2,1,iRI)
    hold on
    plot(orientationLims,[0 0], 'k')
    
    for iCV = 1:3
        plot(groupData.nonTargetOrients{iCV,iRI}, groupData.errors{iCV,iRI}, '.', 'Color', colors{iCV})
    end
    
    set(gca,'XTick', orientationXTicks)
    xlim(orientationLims)
    ylim(errorLims)
    ylabel('error')
    title(targetNames{iRI})
end
xlabel('non-target orientation')
rd_supertitle(sprintf('%s ', subjectIDs{:}));
rd_raiseAxis(gca);

% orientation difference between targets
figure
for iRI = 1:2
    subplot(2,1,iRI)
    hold on
    plot(errorLims,[0 0], 'k')
    
    for iCV = 1:3
        plot(groupData.targetOrientDiff{iCV,iRI}, groupData.errors{iCV,iRI}, '.', 'Color', colors{iCV})
    end
    
    for iCV = 1:3
        smoothError = filter(b,a,errorByODAbs{iCV,iRI});
        plot(odx{iCV,iRI}, smoothError,'Color', colors{iCV}, 'LineWidth',2)
    end
    
    set(gca,'XTick', errorXTicks)
    xlim(errorLims)
    ylim(errorLims)
    ylabel('error')
    title(targetNames{iRI})
end
xlabel('non-target - target orientation difference')
rd_supertitle(sprintf('%s ', subjectIDs{:}));
rd_raiseAxis(gca);

% flip flip
binEdgesFF = -90:10:0;
effVar = [];
figure
for iRI = 1:2
    subplot(2,1,iRI)
    hold on
    for iCV = 1:3
        tod = groupData.targetOrientDiff{iCV,iRI};
        e = groupData.errors{iCV,iRI};
        w = tod>0;
        targetOrientDiffFF{iCV,iRI} = tod;
        targetOrientDiffFF{iCV,iRI}(w) = -tod(w);
        errorsFF{iCV,iRI} = e;
        errorsFF{iCV,iRI}(w) = -e(w);
        
        eff = errorsFF{iCV,iRI}(:);
        todff = targetOrientDiffFF{iCV,iRI}(:);
        
        plot(todff, abs(eff), '.', 'Color', colors{iCV})
        
        % bin
        for iBin = 1:numel(binEdgesFF)-1
            edges = binEdgesFF(iBin:iBin+1);
            effVar(iBin,iCV,iRI) = var(eff(todff>edges(1) & todff<=edges(2)));
        end
    end
    plot([-90 0], [0 0], 'k', 'LineWidth', 2)
    xlim([-90 0])
    xlabel('flipflip non-target - target orientation difference')
    ylabel('flipflip error')
end
todErrorsFF = errorsFF;

figure
for iRI = 1:2
    subplot(1,2,iRI)
    plot(effVar(:,:,iRI))
    xlabel('flipflip non-target - target bin')
    ylabel('error variance')
    title(targetNames{iRI})
end
legend('valid','invalid','neutral')



% target x non-target orientation
validityNames = {'valid','invalid','neutral'};
validityOrder = [1 3 2];
figure
for iRI = 1:2
    for iCV = 1:3
        subplot(3,2,(validityOrder(iCV)-1)*2 + iRI)
        e = errorByTNTO{iCV,iRI};
        e(isnan(e)) = 0;
        sm = smooth2a(e,3,3);
        imagesc(sm')
        colormap(othercolor('PRGn5',256));
        set(gca,'clim',[-8 8])
        axis square
        title(validityNames{iCV})
    end
end

figure
e = allErrorByTNTO;
e(isnan(e)) = 0;
sm = smooth2a(e,3,3);
imagesc(sm')
colormap(othercolor('PRGn5',256));
set(gca,'clim',[-3.5 3.5])
set(gca,'XTickLabel',[])
set(gca,'YTickLabel',[])
axis square
colorbar
xlabel('target orientation')        
ylabel('non-target orientation')        
title('all trials')

if analyzeProbe
    % probe orientation
    figure
    for iRI = 1:2
        subplot(2,1,iRI)
        hold on
        plot(orientationLims,[0 0], 'k')
        
        for iCV = 1:3
            plot(groupData.probeOrients{iCV,iRI}, groupData.errors{iCV,iRI}, '.', 'Color', colors{iCV})
        end
        
        for iCV = 1:3
            smoothError = filter(b,a,errorByPO{iCV,iRI});
            plot(pox{iCV,iRI}, smoothError,'Color', colors{iCV}, 'LineWidth',2)
        end
        
        set(gca,'XTick', orientationXTicks)
        xlim(orientationLims)
        ylim(errorLims)
        ylabel('error')
        title(targetNames{iRI})
    end
    xlabel('probe orientation')
    rd_supertitle(sprintf('%s ', subjectIDs{:}));
    
    % orientation difference between probe and target
    figure
    for iRI = 1:2
        subplot(2,1,iRI)
        hold on
        plot(errorLims,[0 0], 'k')
        
        for iCV = 1:3
            plot(groupData.probeOrientDiff{iCV,iRI}, groupData.errors{iCV,iRI}, '.', 'Color', colors{iCV})
        end
        
        for iCV = 1:3
            smoothError = filter(b,a,errorByPD{iCV,iRI});
            plot(pdx{iCV,iRI}, smoothError,'Color', colors{iCV}, 'LineWidth',2)
        end
        
        set(gca,'XTick', errorXTicks)
        xlim(errorLims)
        ylim(errorLims)
        ylabel('error')
        title(targetNames{iRI})
    end
    xlabel('probe - target orientation difference')
    rd_supertitle(sprintf('%s ', subjectIDs{:}));
    rd_raiseAxis(gca);
end

% smoothed difference between target and non-target orientation
colors = {'b','g','r'};
figure
for iRI = 1:2
    subplot(2,1,iRI)
    hold on
    for iCV = 1:3
        shadedErrorBar(steps, ...
            groupMean.targetOrientDiffSmooth(:,iRI,iCV), ...
            groupSte.targetOrientDiffSmooth(:,iRI,iCV), ...
            colors{iCV},1)
    end
    plot([-90 90],[0 0],'k')
    xlim([steps(1) steps(end)])
    ylim([-15 15])
    ylabel(sprintf('error (sliding window average, size=%d)', winSize))
    if iRI==1
%         legend('valid','invalid','neutral')
    end
end
xlabel('non-target - target orientation difference')
rd_supertitle(sprintf('%s ', subjectIDs{:}));
rd_raiseAxis(gca);

%% Fit lines to data - target orientation
xgridFF = 0:90;
for iSubject = 1:nSubjects
    figure
    hold on
    for iRI = 1:2
        for iCV = 1:3
            eFF = toErrorsFF{iCV,iRI}(:,iSubject);
            toFF = targetOrientsFF{iCV,iRI}(:,iSubject);
            
            p = polyfit(toFF, eFF, 1);
            
            subplot(3,2,(validityOrder(iCV)-1)*2 + iRI)
            hold on
            plot(toFF,eFF,'.')
            plot(xgridFF,polyval(p,xgridFF),'r')
            title(validityNames{iCV})
            
            params{iCV,iRI}(iSubject,:) = p;
        end
    end
    rd_supertitle(subjectIDs{iSubject})
end

for iRI = 1:2
    for iCV = 1:3
        paramsMean(iCV,iRI,:) = mean(params{iCV,iRI},1);
        paramsSte(iCV,iRI,:) = std(params{iCV,iRI},0,1)./sqrt(nSubjects);
    end
end

% plot bars
paramNames = {'slope','intercept'};
ylims = [-.5 .5; -10 25];
for iP = 1:numel(p)
    figure
    for iRI = 1:2
        for iCV = 1:3
            subplot(3,2,(validityOrder(iCV)-1)*2 + iRI)
            bar(params{iCV,iRI}(:,iP))
            title(validityNames{iCV})
            ylim(ylims(iP,:))
        end
    end
    rd_supertitle(paramNames{iP})
end

for iP = 1:numel(p)
    figure
    barweb(paramsMean(validityOrder,:,iP)',paramsSte(validityOrder,:,iP)', ...
        [], targetNames, [], [], [], gray)
    legend(validityNames{validityOrder})
    ylabel(paramNames{iP})
end

fit.to.params = params;
fit.to.paramsMean = paramsMean;
fit.to.paramsSte = paramsSte;
fit.to.paramNames = paramNames;

%% Fit lines to data
xgridFF = -90:0;
for iSubject = 1:nSubjects
    figure
    hold on
    for iRI = 1:2
        for iCV = 1:3
            eFF = todErrorsFF{iCV,iRI}(:,iSubject);
            todFF = targetOrientDiffFF{iCV,iRI}(:,iSubject);
            
            p = polyfit(todFF, eFF, 1);
            
            subplot(3,2,(validityOrder(iCV)-1)*2 + iRI)
            hold on
            plot(todFF,eFF,'.')
            plot(xgridFF,polyval(p,xgridFF),'r')
            title(validityNames{iCV})
            
            params{iCV,iRI}(iSubject,:) = p;
        end
    end
    rd_supertitle(subjectIDs{iSubject})
end

for iRI = 1:2
    for iCV = 1:3
        paramsMean(iCV,iRI,:) = mean(params{iCV,iRI},1);
        paramsSte(iCV,iRI,:) = std(params{iCV,iRI},0,1)./sqrt(nSubjects);
    end
end

% plot bars
paramNames = {'slope','intercept'};
ylims = [-.5 .5; -10 25];
for iP = 1:numel(p)
    figure
    for iRI = 1:2
        for iCV = 1:3
            subplot(3,2,(validityOrder(iCV)-1)*2 + iRI)
            bar(params{iCV,iRI}(:,iP))
            title(validityNames{iCV})
            ylim(ylims(iP,:))
        end
    end
    rd_supertitle(paramNames{iP})
end

for iP = 1:numel(p)
    figure
    barweb(paramsMean(validityOrder,:,iP)',paramsSte(validityOrder,:,iP)', ...
        [], targetNames, [], [], [], gray)
    legend(validityNames{validityOrder})
    ylabel(paramNames{iP})
end

fit.tod.params = params;
fit.tod.paramsMean = paramsMean;
fit.tod.paramsSte = paramsSte;
fit.tod.paramNames = paramNames;
    
%% Set figure properties
% set font size of titles, axis labels, and legends
% set(findall(gcf,'type','text'),'FontSize',14)

