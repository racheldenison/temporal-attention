% rd_plotTemporalAttentionMultiSOAResampledNull.m

%% Setup
subjectInits = {'rd','hl','ho','vp'};
nSubjects = numel(subjectInits);
groupStr = sprintf('N=%d', nSubjects);
nSubjects = numel(subjectInits);

expName = 'E2'; % 'E2','E4'

contrast = 64;

figPrefix = sprintf('g%s_N%d', expName, nSubjects);

saveFigs = 0;

R = load(sprintf('data/E2_SOA_cbD6_randomizationTest_workspace_run98_N%d_20160125.mat',nSubjects));
if ~isequal(subjectInits,R.subjectInits)
    error('different subjects in randomized data')
end

%% Get indiv subject data
for iSubject = 1:nSubjects
    subjectInit = subjectInits{iSubject};
    [acc rt t1t2soa p dp eff] = ...
        rd_plotTemporalAttentionMultiSOA(subjectInit, contrast);
    for iT = 1:numel(acc)
        accData{iT}(:,:,iSubject) = acc{iT};
        rtData{iT}(:,:,iSubject) = rt{iT};
    end
end

%% Group mean
for iT = 1:numel(acc)
    accMean{iT} = mean(acc{iT},3);
    rtMean{iT} = mean(rt{iT},3);
end

%% Pairwise differences
for iT = 1:2
    accDataPairwise(1,iT,:,:) = accData{iT}(1,:,:) - accData{iT}(2,:,:); % VI
    accDataPairwise(2,iT,:,:) = accData{iT}(1,:,:) - accData{iT}(3,:,:); % VN
    accDataPairwise(3,iT,:,:) = accData{iT}(3,:,:) - accData{iT}(2,:,:); % NI
    
    rtDataPairwise(1,iT,:,:) = rtData{iT}(1,:,:) - rtData{iT}(2,:,:);
    rtDataPairwise(2,iT,:,:) = rtData{iT}(1,:,:) - rtData{iT}(3,:,:);
    rtDataPairwise(3,iT,:,:) = rtData{iT}(3,:,:) - rtData{iT}(2,:,:);
end

accMeanPairwise = mean(accDataPairwise,4);
rtMeanPairwise = mean(rtDataPairwise,4);

%% Percentiles for null distributions
p = prctile(R.accDataPairwise,[2.5 97.5],3);
p1 = prctile(R.accDataPairwise,95,3);
pg = prctile(R.accGroupPairwise,[2.5 97.5],3);

%% Stats
[h psig ci stats] = ttest(accDataPairwise,[],[],[],4);

%% Figures
vcNames = {'VI','VN','NI'};
vcColors = {'b','r','r'};
xlims = [0 1000];
for iSubject = 1:nSubjects
    figure
    for iT = 1:2
        for iVC = 1:3
            subplot(3,2,2*(iVC-1)+iT)
            hold on
            plot(xlims,[0 0],':k')
            plot(t1t2soa, squeeze(accDataPairwise(iVC,iT,:,iSubject)),'.-','Color',vcColors{iVC})
            plot(t1t2soa, squeeze(p(iVC,iT,:,:,iSubject))','k')
%             plot(t1t2soa, squeeze(p1(iVC,iT,:,:,iSubject))','color',[.5 .5 .5])
            xlim(xlims)
            if iT==1
                ylabel(vcNames{iVC})
            end
            if iVC==1
                title(sprintf('T%d',iT))
            end
        end
    end
    rd_supertitle(subjectInits{iSubject});
    rd_raiseAxis(gca);
end

% all subjects one one figure, VI only
figure
for iT = 1:2
    for iSubject = 1:nSubjects
        subplot(nSubjects,2,2*(iSubject-1)+iT)
        hold on
        plot(xlims,[0 0],':k')
        plot(t1t2soa, squeeze(accDataPairwise(1,iT,:,iSubject)),'.-','Color',vcColors{1})
        plot(t1t2soa, squeeze(p(1,iT,:,:,iSubject))','k')
%         plot(t1t2soa, squeeze(p1(1,iT,:,:,iSubject))','color',[.5 .5 .5])
        xlim(xlims)
        if iT==1
            ylabel(subjectInits{iSubject})
        end
        if iVC==1
            title(sprintf('T%d',iT))
        end
    end
end
rd_supertitle(vcNames{1});
rd_raiseAxis(gca);

% group data with fixed effects null distribution
ylims = [-0.1 0.1];
figure
for iT = 1:2
    for iVC = 1:3
        subplot(3,2,2*(iVC-1)+iT)
        hold on
        plot(xlims,[0 0],':k')
        plot(t1t2soa, squeeze(accMeanPairwise(iVC,iT,:)),'.-','Color',vcColors{iVC})
        plot(t1t2soa, squeeze(pg(iVC,iT,:,:))','k')
        xlim(xlims)
        ylim(ylims)
        if iT==1
            ylabel(vcNames{iVC})
        end
        if iVC==1
            title(sprintf('T%d',iT))
        end
    end
end
rd_supertitle(sprintf('group, N=%d',nSubjects));
rd_raiseAxis(gca);

figure
for iT = 1:2
    subplot(1,2,iT)
    hold on
    plot(xlims,[0 0],':k')
    plot(t1t2soa, squeeze(accMeanPairwise(1,iT,:)),'.-','Color',vcColors{1})
    plot(t1t2soa, squeeze(pg(1,iT,:,:))','k')
    plot(t1t2soa, squeeze(psig(1,iT,:))','g')
    xlim(xlims)
    ylim(ylims)
    if iT==1
        ylabel(vcNames{1})
    end
    if iVC==1
        title(sprintf('T%d',iT))
    end
end
rd_supertitle(sprintf('group, N=%d',nSubjects));
rd_raiseAxis(gca);

