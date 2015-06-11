% rd_plotTemporalAttentionAdjustBigErrors.m

%% setup
subjectIDs = {'bl','rd','id','ec','ld','en','sj','ml','ca','jl','ew','jx'};
nSubjects = numel(subjectIDs);

run = 9;

errorThresh = 50;

targetNames = {'T1','T2'};

figure
for iSubject = 1:nSubjects
    subjectID = subjectIDs{iSubject};
    [errors, targetOrients, nonTargetOrients, targetOrientDiff, probeOrients, probeOrientDiff, responses, targetOrientDiffSmooth] = ...
        rd_plotTemporalAttentionAdjustErrors(subjectID, run, 0);
    
    ntoErr = [];
    for iRI = 1:2
        subplot(1,2,iRI)
        hold on
        for iCV = 1:3
            nto = nonTargetOrients{iCV,iRI};
            err = errors{iCV,iRI};
            to = targetOrients{iCV,iRI};
            report = to + err;
            
            report(report > 180) = report(report > 180) - 180;
            report(report < 0) = report(report < 0) + 180;
            
            bigErr = abs(err)>errorThresh;
            plot(nto(bigErr), report(bigErr),'.')
            
            ntoErr2 = report(bigErr)-nto(bigErr);
            ntoErr2(ntoErr2>90) = ntoErr2(ntoErr2>90) - 180;
            ntoErr2(ntoErr2<-90) = ntoErr2(ntoErr2<-90) + 180;
            
            ntoError2{iCV,iRI,iSubject} = ntoErr2;
            
            ntoErr = [ntoErr; ntoErr2];
        end
        xlabel('non-target orientation')
        if iRI==1
            ylabel('reported orientation')
        end
        title(targetNames{iRI})
    end
    ntoError{iSubject} = ntoErr;
end

%% aggregate and plot data
binSize = 10;
x = -90:binSize:90;

% all conditions combined
figure
iS = 0;
for i = 1:4
    for j = 1:3
        iS = iS + 1;
        subplot(4,3,iS)
        count(:,iS) = hist(ntoError{iS}, x);
        bar(x,count(:,iS));
    end
end
countMean = mean(count,2);
countSte = std(count,0,2)./sqrt(nSubjects);

figure
shadedErrorBar(x, countMean, countSte);
xlim([-90 90])
xlabel('difference from non-target orientation')
ylabel('number of reports')

% by condition
for iS = 1:nSubjects
    for iRI = 1:2
        for iCV = 1:3
            countByCond(:,iCV,iRI,iS) = hist(ntoError2{iCV,iRI,iS}, x);
        end
    end
end
countByCondMean = mean(countByCond,4);

for iS = 1:nSubjects
    for iRI = 1:2
        for iCV = 1:3
            rateByCond(:,iCV,iRI,iS) = countByCond(:,iCV,iRI,iS)./length(targetOrients{iCV,iRI});
        end
    end
end
rateByCondMean = mean(rateByCond,4);
rateByCondSte = std(rateByCond,0,4)./sqrt(nSubjects);

colors = {'b','g','r'};
figure
for iRI = 1:2
    subplot(1,2,iRI)
    hold on
    for iCV = 1:3
        %     plot(x, rateByCondMean(:,:,iRI))
        shadedErrorBar(x, rateByCondMean(:,iCV,iRI), rateByCondSte(:,iCV,iRI), colors{iCV}, 1);
    end
    title(targetNames{iRI})
    xlim([-90 90])
    ylim([0 12e-3])
    set(gca,'XTick',[-90 0 90])
end
