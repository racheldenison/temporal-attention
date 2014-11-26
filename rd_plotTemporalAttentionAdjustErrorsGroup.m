% rd_plotTemporalAttentionAdjustErrorsGroup.m

%% setup
subjectIDs = {'bl','rd','id','ec','ld'};
nSubjects = numel(subjectIDs);

analyzeProbe = 0;

%% get data
for iSubject = 1:nSubjects
    subjectID = subjectIDs{iSubject};
    [groupData0(iSubject).errors, ...
        groupData0(iSubject).targetOrients, ...
        groupData0(iSubject).nonTargetOrients, ...
        groupData0(iSubject).targetOrientDiff, ...
        groupData0(iSubject).probeOrients, ...
        groupData0(iSubject).probeOrientDiff] = ...
            rd_plotTemporalAttentionAdjustErrors(subjectID);
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
                groupData.(fName){iCV,iRI}(:,iSubject) = ...
                    groupData0(iSubject).(fName){iCV,iRI};
            end
        end
    end
end

% calculate mean errors for each x-axis value
for iRI = 1:2
    for iCV = 1:3
        allErrors = groupData.errors{iCV,iRI}(:);
        allTargetOrients = groupData.targetOrients{iCV,iRI}(:);
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
        pox{iCV,iRI} = pos;
        pdx{iCV,iRI} = pds;
    end
end

%% plot figures
targetNames = {'T1','T2'};
colors = {'b','g','r'};
errorLims = [-100 100];
orientationLims = [-10 190];
errorXTicks = [-90 -45 0 45 90];
orientationXTicks = [0 45 90 135 180];

smoothSize = 5;
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
