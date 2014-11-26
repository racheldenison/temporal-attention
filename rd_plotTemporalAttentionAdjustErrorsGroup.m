% rd_plotTemporalAttentionAdjustErrorsGroup.m

%% setup
subjectIDs = {'bl','rd','id','ec'};
nSubjects = numel(subjectIDs);

%% get data
for iSubject = 1:nSubjects
    subjectID = subjectIDs{iSubject};
    [groupData0(iSubject).errors, ...
        groupData0(iSubject).targetOrients, ...
        groupData0(iSubject).nonTargetOrients, ...
        groupData0(iSubject).targetOrientDiff] = ...
            rd_plotTemporalAttentionAdjustErrors(subjectID);
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

        ods = unique(allTargetOrientDiff);
        for iOD = 1:numel(ods)
            odIdx = allTargetOrientDiff==ods(iOD);
            errorByOD{iCV,iRI}(iOD) = mean(allErrors(odIdx));
        end
        
        tos = unique(allTargetOrients);
        for iTO = 1:numel(tos)
            toIdx = allTargetOrients==tos(iTO);
            errorByTO{iCV,iRI}(iTO) = mean(allErrors(toIdx));
        end
        
        odx{iCV,iRI} = ods;
        tox{iCV,iRI} = tos;
    end
end

%% plot figures
targetNames = {'T1','T2'};
colors = {'b','g','r'};
errorLims = [-100 100];
orientationLims = [-10 190];
errorXTicks = [-90 -45 0 45 90];
orientationXTicks = [0 45 90 135 180];

smoothSize = 3;
b = (1/smoothSize)*ones(1,smoothSize);
a = 1;

% target orientation
figure
for iRI = 1:2
    subplot(2,1,iRI)
    hold on
    plot(orientationLims,[0 0], 'k')
    
    for iCV = 1:3
        plot(groupData.targetOrients{iCV,iRI}, groupData.errors{iCV,iRI}, '.', 'Color', colors{iCV})
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
        smoothError = filter(b,a,errorByOD{iCV,iRI});
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
