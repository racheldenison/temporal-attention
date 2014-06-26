function [accMean rtMean t1t2soa p dprime crit eff] = rd_plotTemporalAttentionMultiSOA(subjectInit)

subjectInit = 'rd';
exptName = 'cbD6';
% tilt = '*';
contrast = 64;

soa1 = [1000 1000 1000 1000];
soa2 = [1100 1150 1200 1250];
t1t2soa = soa2 - soa1;
run = 9;

cleanRT = 0;

dataDir = pathToExpt('data');

% subjectID = sprintf('%s_%s_tilt%s_tc%d', ...
%     subjectInit, exptName, tilt, contrast);
subjectID = sprintf('%s_%s*tc%d', ...
    subjectInit, exptName, contrast);

if cleanRT
    rtStr = '_RTx';
else
    rtStr = '';
end

%% Get data
for iSOA = 1:numel(soa1)
    subject = sprintf('%s_soa%d-%d', ...
        subjectID, soa1(iSOA), soa2(iSOA));
    
    % load data from a given soa
    dataFile = dir(sprintf('%s/%s_run%02d%s_T*', ...
        dataDir, subject, run, rtStr));
    if numel(dataFile)~=1
        fprintf('%s/%s_run%02d%s_T*', dataDir, subject, run, rtStr)
        error('more or fewer than one matching data file')
    else
        load(sprintf('%s/%s', dataDir, dataFile.name))
    end
    
    % read out the target timings
    soas(iSOA,:) = expt.p.soas;
    
    % check soas
    if soas(iSOA,1)~=soa1(iSOA)/1000 || soas(iSOA,2)~=soa2(iSOA)/1000
        error('requested SOAs do not match SOAs in data file')
    end
    
    % read out the accuracy and rt
    for iEL = 1:2 % early/late
        accMean{iEL}(:,iSOA) = results.accMean{iEL};
        accSte{iEL}(:,iSOA) = results.accSte{iEL};
        rtMean{iEL}(:,iSOA) = results.rtMean{iEL};
        rtSte{iEL}(:,iSOA) = results.rtSte{iEL};
    end
    
    % calculate dprime
    results = rd_analyzeTemporalAttentionDprime(expt, results);
    % read out dprime and crit
    for iEL = 1:2 % early/late
        dprime{iEL}(:,iSOA) = results.totals.dprime{iEL};
        crit{iEL}(:,iSOA) = results.totals.crit{iEL};
    end
end

p = expt.p;

%% Calculate efficiency
for iEL = 1:2 % early/late
    eff{iEL} = dprime{iEL}./rtMean{iEL};
end

%% Plot figs
intervalNames = {'early','late'};
accLims = [0.2 1];
rtLims = [0.2 1.6];
dpLims = [-0.5 2.7];
critLims = [-1 1];
effLims = [-0.5 8];
soaLims = [t1t2soa(1)-100 t1t2soa(end)+100];
colors = get(0,'DefaultAxesColorOrder');
axTitle = '';

fig(1) = figure;
for iRI = 1:numel(p.respInterval)
    subplot(1,numel(p.respInterval),iRI)
    hold on
    plot(soaLims, [0.5 0.5], '--k');
    
    p1 = errorbar(repmat(t1t2soa',1,numel(p.cueValidity)),...
        accMean{iRI}', accSte{iRI}', '.-', 'MarkerSize', 20);

    xlabel('soa')
    ylabel('acc')
    legend(p1, num2str(p.cueValidity'),'location','best')
    title(intervalNames{iRI})
    xlim(soaLims)
    ylim(accLims)
    rd_supertitle(subjectID);
    rd_raiseAxis(gca);
    rd_supertitle(axTitle);
end

fig(2) = figure;
for iRI = 1:numel(p.respInterval)
    subplot(1,numel(p.respInterval),iRI)
    
    p1 = errorbar(repmat(t1t2soa',1,numel(p.cueValidity)),...
        rtMean{iRI}', rtSte{iRI}', '.-', 'MarkerSize', 20);
    
    xlabel('soa')
    ylabel('rt')
    legend(num2str(p.cueValidity'),'location','best')
    title(intervalNames{iRI})
    xlim(soaLims)
    ylim(rtLims)
    box off
    rd_supertitle(subjectID);
    rd_raiseAxis(gca);
    rd_supertitle(axTitle);
end

fig(3) = figure;
for iRI = 1:numel(p.respInterval)
    subplot(1,numel(p.respInterval),iRI)
    hold on
    plot(soaLims, [0 0], '--k');
    
    p1 = plot(repmat(t1t2soa',1,numel(p.cueValidity)),...
        dprime{iRI}', '.-', 'MarkerSize', 20);

    xlabel('soa')
    ylabel('dprime')
    legend(p1, num2str(p.cueValidity'),'location','best')
    title(intervalNames{iRI})
    xlim(soaLims)
    ylim(dpLims)
    rd_supertitle(subjectID);
    rd_raiseAxis(gca);
    rd_supertitle(axTitle);
end

fig(4) = figure;
for iRI = 1:numel(p.respInterval)
    subplot(1,numel(p.respInterval),iRI)
    hold on
    plot(soaLims, [0 0], '--k');
    
    p1 = plot(repmat(t1t2soa',1,numel(p.cueValidity)),...
        crit{iRI}', '.-', 'MarkerSize', 20);

    xlabel('soa')
    ylabel('crit')
    legend(p1, num2str(p.cueValidity'),'location','best')
    title(intervalNames{iRI})
    xlim(soaLims)
    ylim(critLims)
    rd_supertitle(subjectID);
    rd_raiseAxis(gca);
    rd_supertitle(axTitle);
end

fig(5) = figure;
for iRI = 1:numel(p.respInterval)
    subplot(1,numel(p.respInterval),iRI)
    hold on
    plot(soaLims, [0 0], '--k');
    
    p1 = plot(repmat(t1t2soa',1,numel(p.cueValidity)),...
        eff{iRI}', '.-', 'MarkerSize', 20);

    xlabel('soa')
    ylabel('eff')
    legend(p1, num2str(p.cueValidity'),'location','best')
    title(intervalNames{iRI})
    xlim(soaLims)
    ylim(effLims)
    rd_supertitle(subjectID);
    rd_raiseAxis(gca);
    rd_supertitle(axTitle);
end
