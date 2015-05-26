% rd_fitDoG

%% setup
subjectIDs = {'bl','rd','id','ec','ld','en','sj','ml','ca','jl','ew','jx'};
nSubjects = numel(subjectIDs);

run = 9;
plotFigs = 1;
plotIndivFigs = 0;

p0 = [50 0];

xgrid = -90:90;

%% get data and fit
for iSubject = 1:nSubjects
    subjectID = subjectIDs{iSubject};
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

%% plot figs
if plotFigs
    for iSubject = 1:nSubjects
        figure
        for iRI = 1:2
            for iCV = 1:3
                subplot(3,2,(iCV-1)*2 + iRI)
                hold on
                p = params{iCV,iRI}(iSubject,:);
                plot(groupData(iSubject).tod{iCV,iRI}, groupData(iSubject).e{iCV,iRI},'.')
                plot(xgrid,DoG(0,p(1),p(2),0,xgrid),'r')
            end
        end
        rd_supertitle(subjectIDs{iSubject});
    end
end

paramNames = {'sd','amp'};
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
        
validityOrder = [1 3 2];
figure
barweb(paramsMean(validityOrder,:,2)',paramsSte(validityOrder,:,2)', ...
    [], {'T1','T2'}, [], [], [], gray)
legend('valid','neutral','invalid')
ylabel('DoG amplitude')



