% rd_fitDoG

%% setup
subjectIDs = {'bl','rd','id','ec','ld','en','sj','ml','ca','jl','ew','jx'};
nSubjects = numel(subjectIDs);

run = 9;
plotFigs = 1;
plotIndivFigs = 0;

% p0 = [0 50];
p0 = 0;

xgrid = -90:90;

%% get data and fit
for iSubject = 1:nSubjects
    subjectID = subjectIDs{iSubject};
    fprintf('%s\n', subjectID)
    [e, to, nto, tod] = ...
        rd_plotTemporalAttentionAdjustErrors(subjectID, run, plotIndivFigs);

    for iRI = 1:2
        for iCV = 1:3
            xdata = tod{iCV,iRI};
            ydata = e{iCV,iRI};
            
            p = fminsearch(@(p) rd_DoGCost(p,xdata,ydata), p0);
            
            params{iCV,iRI}(iSubject,:) = p;
        end
    end
    
    groupData(iSubject).e = e;
    groupData(iSubject).to = to;
    groupData(iSubject).nto = nto;
    groupData(iSubject).tod = tod;
end

for iRI = 1:2
    for iCV = 1:3
        paramsMean(iCV,iRI,:) = mean(params{iCV,iRI},1);
        paramsSte(iCV,iRI,:) = std(params{iCV,iRI},0,1)./sqrt(nSubjects);
    end
end

%% fit the data points from all subjects together
% organize data
for iRI = 1:2
    for iCV = 1:3
        groupDataAll.e{iCV,iRI} = [];
        groupDataAll.tod{iCV,iRI} = [];
        for iSubject = 1:nSubjects
            groupDataAll.e{iCV,iRI} = [groupDataAll.e{iCV,iRI} groupData(iSubject).e{iCV,iRI}];
            groupDataAll.tod{iCV,iRI} = [groupDataAll.tod{iCV,iRI} groupData(iSubject).tod{iCV,iRI}];
        end
    end
end

% fit
for iRI = 1:2
    for iCV = 1:3
        xdata = groupDataAll.tod{iCV,iRI}(:);
        ydata = groupDataAll.e{iCV,iRI}(:);
        
        p = fminsearch(@(p) rd_DoGCost(p,xdata,ydata), p0);
        
        groupParams{iCV,iRI} = p;
    end
end

%% plot figs
targetNames = {'T1','T2'};
validityOrder = [1 3 2];
validityNames = {'valid','invalid','neutral'};

if plotFigs
    for iSubject = 1:nSubjects
        figure
        for iRI = 1:2
            for iCV = 1:3
                subplot(3,2,(iCV-1)*2 + iRI)
                hold on
                p = params{iCV,iRI}(iSubject,:);
                plot(groupData(iSubject).tod{iCV,iRI}, groupData(iSubject).e{iCV,iRI},'.')
                switch numel(p)
                    case 1
                        plot(xgrid,DoG(0,44,p(1),0,xgrid),'r')
                    case 2
                        plot(xgrid,DoG(0,p(2),p(1),0,xgrid),'r')
                    otherwise
                        error('p has wrong number of elements')
                end
            end
        end
        rd_supertitle(subjectIDs{iSubject});
    end
end

switch numel(p)
    case 1
        paramNames = {'amp'};
    case 2
        paramNames = {'amp','sd'};
end
for iP = 1:length(p)
    figure
    for iRI = 1:2
        for iCV = 1:3
            subplot(3,2,(iCV-1)*2 + iRI)
            bar(params{iCV,iRI}(:,iP))
            if iP==2
                ylim([-20 20])
            end
        end
    end
    rd_supertitle(paramNames{iP})
end
        
figure
barweb(paramsMean(validityOrder,:,1)',paramsSte(validityOrder,:,1)', ...
    [], targetNames, [], [], [], gray)
legend(validityNames{validityOrder})
ylabel('DoG amplitude')

% group fit
figure
for iRI = 1:2
    for iCV = 1:3
        subplot(3,2,(validityOrder(iCV)-1)*2 + iRI)
        hold on
        p = groupParams{iCV,iRI};
        plot(groupDataAll.tod{iCV,iRI}(:), groupDataAll.e{iCV,iRI}(:),'.')
        switch numel(p)
            case 1
                plot(xgrid,DoG(0,44,p(1),0,xgrid),'r','LineWidth',2)
            case 2
                plot(xgrid,DoG(0,p(2),p(1),0,xgrid),'r','LineWidth',2)
            otherwise
                error('p has wrong number of elements')
        end
        title(validityNames{iCV})
        if validityOrder(iCV)==3 && iRI==1
            xlabel('non-target - target orientation (deg)')
            ylabel('error (deg)')
        end
    end
end
h = rd_supertitle('group combined data');
% rd_raiseAxis(h);

% organize again
for iRI = 1:2
    for iCV = 1:3
        for iP = 1:numel(p)
            gp.(paramNames{iP})(iCV,iRI) = groupParams{iCV,iRI}(:,iP);
        end
    end
end

figure
for iP = 1:numel(p)
    subplot(1,numel(p),iP)
    bar(gp.(paramNames{iP})(validityOrder,:)')
    title(paramNames{iP})
    set(gca,'XTickLabel',targetNames)
end
colormap gray
legend(validityNames{validityOrder})




