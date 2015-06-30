function [accMean rtMean t1t2soa p dprime crit eff] = rd_plotTemporalAttentionMultiSOA(subjectInit)

% subjectInit = 'vp';
exptName = 'cbD10'; % 'cbD6'
tilt = '*';
contrast = 64; % plot one contrast at a time

% soa1 = [1000 1000 1000 1000 1000 1000 1000 1000 1000 1000];
% soa2 = [1100 1150 1200 1250 1300 1350 1400 1450 1500 1800];
% soa1 = [1000 1000 1000 1000 1000 1000 1000 1000 1000];
% soa2 = [1100 1150 1200 1250 1300 1400 1450 1500 1800];
soa1 = [1000 1000 1000];
soa2 = [1100 1300 1800];
t1t2soa = soa2 - soa1;
run = 8;

plotFigs = 0;

cleanRT = 0;
doDprime = 0;

expName = 'E4_contrast_cbD10'; % 'E2_SOA_cbD6'
dataDir = pathToExpt('data');
dataDir = sprintf('%s/%s/%s', dataDir, expName, subjectInit(1:2));

% subjectID = sprintf('%s_%s_tilt%s_tc%d', ...
%     subjectInit, exptName, tilt, contrast);
subjectID = sprintf('%s_%s*', ...
    subjectInit, exptName);

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
    
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%% if you want to reanalyze, do it here %%%
%     T1T2Axis = 'diff';
%     [expt results] = rd_analyzeTemporalAttention(expt, 0, 0, 0, 0, T1T2Axis, 0);
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % select target contrast
    tcIdx = find(expt.p.targetContrasts==contrast/100);
    if isempty(tcIdx)
        error('requested target contrast not found')
    end

    % read out the target timings
    soas(iSOA,:) = expt.p.soas;
    
    % check soas
    if soas(iSOA,1)~=soa1(iSOA)/1000 || soas(iSOA,2)~=soa2(iSOA)/1000
        error('requested SOAs do not match SOAs in data file')
    end
    
    % read out the accuracy and rt
    for iEL = 1:2 % early/late
        accMean{iEL}(:,iSOA) = results.accMean{iEL}(:,tcIdx);
        accSte{iEL}(:,iSOA) = results.accSte{iEL}(:,tcIdx);
        rtMean{iEL}(:,iSOA) = results.rtMean{iEL}(:,tcIdx);
        rtSte{iEL}(:,iSOA) = results.rtSte{iEL}(:,tcIdx);
        
        if doDprime
            [dprime{iEL}(:,iSOA) crit{iEL}(:,iSOA)] = ...
                rd_dprime(results.accMean{iEL}(:,tcIdx),[],'2afc');
        else
            dprime{iEL}(:,iSOA) = zeros(size(rtSte{iEL}(:,iSOA)));
            crit{iEL}(:,iSOA) = zeros(size(rtSte{iEL}(:,iSOA)));
        end
    end
    
    %     % calculate dprime
    %     results = rd_analyzeTemporalAttentionDprime(expt, results);
    %     % read out dprime and crit
    %     for iEL = 1:2 % early/late
    %         dprime{iEL}(:,iSOA) = results.totals.dprime{iEL};
    %         crit{iEL}(:,iSOA) = results.totals.crit{iEL};
    %     end
end

p = expt.p;

%% Calculate efficiency
for iEL = 1:2 % early/late
    if doDprime
        eff{iEL} = dprime{iEL}./rtMean{iEL};
    else
        eff{iEL}(:,iSOA) = zeros(size(rtSte{iEL}(:,iSOA)));
    end
end

%% Plot figs
if plotFigs
    intervalNames = {'early','late'};
    accLims = [0.2 1];
    switch exptName
        case 'cbD10'
            rtLims = [0.8 2.2];
        otherwise
            rtLims = [0.2 1.6];
    end
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
end
