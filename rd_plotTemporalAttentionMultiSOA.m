function [accMean rtMean t1t2soa p dprime eff] = rd_plotTemporalAttentionMultiSOA(subjectInit, contrast, T1T2Axis, extraSelection, cleanRT)

if nargin < 3
    T1T2Axis = []; % 'same','diff'
end
if nargin < 4
    extraSelection = []; % e.g. 'sameOrient','diffOrient'
end
if nargin < 5
    cleanRT = 0;
end

% subjectInit = 'vp';
exptName = 'cbD6'; % 'cbD6', 'cbD10'
tilt = '*';
% contrast = 16; % plot one contrast at a time

soa1 = [1000 1000 1000 1000 1000 1000 1000 1000 1000 1000];
soa2 = [1100 1150 1200 1250 1300 1350 1400 1450 1500 1800];
t1t2soa = soa2 - soa1;
run = 98; % 8 = runs 1-3; 9 = runs 2-3; 18 = runs 1-3 with first good block of each day excluded; 98 = all runs

plotFigs = 1;
saveFigs = 1;

doDprime = 1;

expName = 'E2_SOA_cbD6'; % 'E2_SOA_cbD6', 'E4_contrast_cbD10'
dataDir = pathToExpt('data');
dataDir = sprintf('%s/%s/%s', dataDir, expName, subjectInit(1:2));
figDir = sprintf('%s/%s/%s', pathToExpt('figures'), expName, subjectInit(1:2));

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
    dataFile = dir(sprintf('%s/%s_run%02d%s_T*.mat', ...
        dataDir, subject, run, rtStr));
    if numel(dataFile)~=1
        fprintf('%s/%s_run%02d%s_T*', dataDir, subject, run, rtStr)
        error('more or fewer than one matching data file')
    else
        load(sprintf('%s/%s', dataDir, dataFile.name))
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% if you want to reanalyze, do it here %%%
    if ~isempty(T1T2Axis) || ~isempty(extraSelection)
        [expt results] = rd_analyzeTemporalAttention(expt, 0, 0, 0, 0, T1T2Axis, 0, extraSelection);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
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
            dprime{iEL}(:,iSOA) = ...
                rd_dprime(results.accMean{iEL}(:,tcIdx),[],'2afc','adjust');
        else
            dprime{iEL}(:,iSOA) = zeros(size(rtSte{iEL}(:,iSOA)));
        end
    end
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
    effLims = [-0.5 4];
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
        title(intervalNames{iRI})
        xlim(soaLims)
        ylim(accLims)
        legend(p1, num2str(p.cueValidity'),'location','best')
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
        title(intervalNames{iRI})
        xlim(soaLims)
        ylim(rtLims)
        legend(num2str(p.cueValidity'),'location','best')
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
        title(intervalNames{iRI})
        xlim(soaLims)
        ylim(dpLims)
        legend(p1, num2str(p.cueValidity'),'location','best')
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
            eff{iRI}', '.-', 'MarkerSize', 20);
        
        xlabel('soa')
        ylabel('eff')
        title(intervalNames{iRI})
        xlim(soaLims)
        ylim(effLims)
        legend(p1, num2str(p.cueValidity'),'location','best')
        rd_supertitle(subjectID);
        rd_raiseAxis(gca);
        rd_supertitle(axTitle);
    end
end

%% Save figs
if saveFigs
    if isempty(T1T2Axis)
        T1T2Str = 'T1T2all';
    end
    figPrefix = sprintf('%s_%s_contrast%d_run%02d%s_%s', ...
        subjectInit, exptName, contrast, run, rtStr, T1T2Str);
    figNames = {'acc','rt','dprime','eff'};
    rd_saveAllFigs(fig, figNames, figPrefix, figDir)
end
