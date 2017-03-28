function [pairedBC, pairNames, stackedBC, stackedNames] = rd_plotTemporalAttentionTradeoffsMean2(sample)

%% setup
e0 = load('data/E0_workspace_run09_N10_20160224.mat');
e3 = load('data/E3_workspace_run09_N12_20160224.mat');
e5 = load('data/E5_workspace_run01_N10_20160806.mat');

% fractional benefits and costs?
normalizeWithinTarget = 1;
baselineOpt = 'neut'; % 'vi','condsum','neut'

% for generating pairedBC bootstrap error bars or null distributions with
% rd_resampleTemporalAttentionTradeoffs
% to plot real data, turn these OFF
resampleSubjects = 0; % draws one resample from subject data
resampleTrials = 0; % loads one resample from pre-saved data file
resampleOption = 'bootstrap'; % 'bootstrap','permutation'

% for plotting bootstrap error bars that have been generated previously
% must load below
bootstrapErrorBars = 0;

% for randomization stats
doRandomizationStats = 0;

plotFigs = 0;

%% resample if requested
% currently does not include condsum
if resampleSubjects
    e0.accDataCP0 = e0.accDataCP;
    for iT = 1:size(e0.accDataCP,3)
        for iP = 1:size(e0.accDataCP,1)
            vals = e0.accDataCP(iP,:,iT);
            e0.accDataCP(iP,:,iT) = RandSample(vals,size(vals));
        end
    end
    e3.pdData0 = e3.pdData;
    for iT = 1:size(e3.pdData,3) 
        for iP = 1:size(e3.pdData,1)
            vals = e3.pdData(iP,:,iT);
            e3.pdData(iP,:,iT) = RandSample(vals,size(vals));
        end
    end
    e5.accDataCP0 = e5.accDataCP;
    for iT = 1:size(e5.accDataCP,3)
        for iP = 1:size(e5.accDataCP,1)
            vals = e5.accDataCP(iP,:,iT);
            e5.accDataCP(iP,:,iT) = RandSample(vals,size(vals));
        end
    end
    e5.accDataIBP0 = e5.accDataIBP;
    for iT = 1:size(e5.accDataIBP,3)
        for iP = 1:size(e5.accDataIBP,1)
            vals = e5.accDataIBP(iP,:,iT);
            e5.accDataIBP(iP,:,iT) = RandSample(vals,size(vals));
        end
    end
end

if resampleTrials
    switch resampleOption
        case 'bootstrap'
            % load bootstrap workspace, overwrites previous e0, e3, e5
            e0 = load('data/E0_resample_workspace_bootstrap_20161128.mat');
            accDataCP = e0.accDataCPairwise(:,:,sample);
            e0.accDataCP = reshape(accDataCP,[3,1,2]);
            
            e3 = load('data/E3_resample_workspace_bootstrap_20161128.mat');
            pdData = e3.pd(:,:,2,sample);
            e3.pdData = reshape(pdData,[3,1,2]);
            
            e5 = load('data/E5_resample_workspace_bootstrap_20161128.mat');
            accDataCP = e5.accDataCPairwise(:,:,sample);
            accDataIBP = e5.accDataIBPairwise(:,:,sample);
            e5.accDataCP = reshape(accDataCP,[3,1,3]);
            e5.accDataIBP = reshape(accDataIBP,[5,1,3]);
        case 'permutation'
            % load randomization workspace, overwrites previous e0, e3, e5
            e0 = load('data/E0_randomizationTest_workspace_run09_N10_20160107b.mat');
            accDataCP = e0.accDataCPairwise(:,:,sample);
            e0.accDataCP = reshape(accDataCP,[3,1,2]);
            
            e3 = load('data/adjust_randomizationTest_workspace_run09_N12_20160108.mat');
            pdData = e3.pd(:,:,2,sample);
            e3.pdData = reshape(pdData,[3,1,2]);
            
            e5 = load('data/E5_randomizationTest_workspace_run01_IB_N10_20160822.mat');
            accDataCP = e5.accDataCPairwise(:,:,sample);
            accDataIBP = e5.accDataIBPairwise(:,:,sample);
            e5.accDataCP = reshape(accDataCP,[3,1,3]);
            e5.accDataIBP = reshape(accDataIBP,[5,1,3]);
        otherwise
            error('resampleOption not recognized')
    end
end

%% get data
% valid vs. invalid
vi.e0 = squeeze(mean(e0.accDataCP(1,:,:),2));
vi.e3 = squeeze(mean(e3.pdData(1,:,:),2));
vi.e5 = squeeze(mean(e5.accDataCP(1,:,:),2));

vi12.e5 = squeeze(mean(e5.accDataIBP(1:2,:,:),2))'; % target x VI1/VI2

if ~resampleSubjects && ~resampleTrials
    % sum across conditions
    condsum.e0 = sum(e0.accMeanC)';
    condsum.e3 = sum(e3.paramsMean.sd)';
    condsum.e5 = sum(e5.accMeanC)';
    
    condsum12.e5(:,1) = sum(e5.accMeanIB(1:3,:)); % I1
    condsum12.e5(:,2) = sum(e5.accMeanIB([1 2 4],:)); % I2
    
    % neutral
    neut.e0 = e0.accMeanC(3,:)';
    neut.e3 = e3.paramsMean.sd(3,:)';
    neut.e5 = e5.accMeanC(3,:)';
    
    neut12.e5(:,1) = e5.accMeanIB(2,:)'; % I1, same as neut.e5
    neut12.e5(:,2) = e5.accMeanIB(2,:)'; % I2, same as neut.e5    
end

% set baseline
base = eval(baselineOpt);
base12 = eval([baselineOpt '12']);

% mean benefits and costs
% VI, VN, NI
bc.e0.t1b = mean(e0.accDataCP(2,:,1),2);
bc.e0.t1c_cuet2 = mean(e0.accDataCP(3,:,1),2);
bc.e0.t2b = mean(e0.accDataCP(2,:,2),2);
bc.e0.t2c_cuet1 = mean(e0.accDataCP(3,:,2),2);

bc.e3.t1b = mean(e3.pdData(2,:,1),2);
bc.e3.t1c_cuet2 = mean(e3.pdData(3,:,1),2);
bc.e3.t2b = mean(e3.pdData(2,:,2),2);
bc.e3.t2c_cuet1 = mean(e3.pdData(3,:,2),2);

% VI1, VI2, VN, NI1, NI2
bc.e5.t1b = mean(e5.accDataIBP(3,:,1),2);
bc.e5.t2b = mean(e5.accDataIBP(3,:,2),2);
bc.e5.t3b = mean(e5.accDataIBP(3,:,3),2);
bc.e5.t1c_cuet2 = mean(e5.accDataIBP(4,:,1),2);
bc.e5.t1c_cuet3 = mean(e5.accDataIBP(5,:,1),2);
bc.e5.t2c_cuet1 = mean(e5.accDataIBP(4,:,2),2);
bc.e5.t2c_cuet3 = mean(e5.accDataIBP(5,:,2),2);
bc.e5.t3c_cuet1 = mean(e5.accDataIBP(4,:,3),2);
bc.e5.t3c_cuet2 = mean(e5.accDataIBP(5,:,3),2);

% % benefits and costs: effect size (Cohen's d)
% % VI, VN, NI
% bc.e0.t1b = e0.dP(2,:,1);
% bc.e0.t1c_cuet2 = e0.dP(3,:,1);
% bc.e0.t2b = e0.dP(2,:,2);
% bc.e0.t2c_cuet1 = e0.dP(3,:,2);
% 
% bc.e3.t1b = e3.dP(2,:,1);
% bc.e3.t1c_cuet2 = e3.dP(3,:,1);
% bc.e3.t2b = e3.dP(2,:,2);
% bc.e3.t2c_cuet1 = e3.dP(3,:,2);
% 
% % VI1, VI2, VN, NI1, NI2
% bc.e5.t1b = e5.dP(3,:,1);
% bc.e5.t2b = e5.dP(3,:,2);
% bc.e5.t3b = e5.dP(3,:,3);
% bc.e5.t1c_cuet2 = e5.dP(4,:,1);
% bc.e5.t1c_cuet3 = e5.dP(5,:,1);
% bc.e5.t2c_cuet1 = e5.dP(4,:,2);
% bc.e5.t2c_cuet3 = e5.dP(5,:,2);
% bc.e5.t3c_cuet1 = e5.dP(4,:,3);
% bc.e5.t3c_cuet2 = e5.dP(5,:,3);

expNames = {'e0','e3','e5'};
nExp = numel(expNames);

targetNames.e0 = {'t1','t2'};
targetNames.e3 = {'t1','t2'};
targetNames.e5 = {'t1','t2','t3'};

%% pair targets
for iExp = 1:nExp
    expName = expNames{iExp};
    tn = targetNames.(expName);
    alltargets = 1:numel(tn);
    for iT = 1:numel(tn)
        nontargets = setdiff(alltargets,iT);
        for iNT = 1:numel(nontargets)
            bn = sprintf('%sb', tn{iT});
            benefit = bc.(expName).(bn);
            nt = nontargets(iNT);
            cn = sprintf('%sc_cue%s',tn{nt},tn{iT});
            cost = bc.(expName).(cn);
            if normalizeWithinTarget
                if strcmp(expName,'e5')
                    targets = setdiff(alltargets,nt);
                    idx = find(targets==iT);
                    benefit = benefit/base12.(expName)(iT,iNT);
                    cost = cost/base12.(expName)(nt,idx);
                else
                    benefit = benefit/base.(expName)(iT);
                    cost = cost/base.(expName)(nt);
                end
            end
            pairedBC.(expName)(iT,:,iNT) = [benefit cost];
            pairNames.(expName){iT,iNT} = sprintf('%s_%s',bn,cn);
        end
    end
end

%% stack targets
for iExp = 1:nExp
    expName = expNames{iExp};
    tn = targetNames.(expName);
    alltargets = 1:numel(tn);
    for iT = 1:numel(tn)
        nontargets = setdiff(alltargets,iT);
        for iNT = 1:numel(nontargets)
            bn = sprintf('%sb', tn{iT});
            benefit = bc.(expName).(bn);
            nt = nontargets(iNT);
            cn = sprintf('%sc_cue%s',tn{iT},tn{nt});
            cost = bc.(expName).(cn);
            if normalizeWithinTarget
                if strcmp(expName,'e5')
                    benefit = benefit/base12.(expName)(iT,iNT);
                    cost = cost/base12.(expName)(iT,iNT);
                else
                    benefit = benefit/base.(expName)(iT);
                    cost = cost/base.(expName)(iT);
                end
            end
            stackedBC.(expName)(iT,:,nt) = [benefit cost];
            stackedNames.(expName){iT,nt} = sprintf('%s_%s',bn,cn);
        end
    end
end

%% plot figs
if plotFigs
%% plot scatterplots
colors = {'b','r','k'};
% colors = {[97 47 255]/255,[255 4 0]/255,[255 177 6]/255};
% colors = {[.1 .1 .1], [.4 .4 .4], [.7 .7 .7]};
% shapes = {'o','s','^'};
shapes = {'o','^','s'};
faceColors = {[1 1 1],[.5 .5 .5]};
faceColors2 = colors;
figure
hold on
for iExp = 1:nExp
    expName = expNames{iExp};
    dataMax = max(pairedBC.(expName)(:));
    dataMean = mean(pairedBC.(expName)(:));
    dataSD = std(pairedBC.(expName)(:));
    for iNT = 1:size(pairedBC.(expName),3)
        vals = pairedBC.(expName)(:,:,iNT);
%         vals = pairedBC.(expName)(:,:,iNT)/dataMax;
%         vals = (pairedBC.(expName)(:,:,iNT)-dataMean)/dataSD;
        plot(vals(:,1),vals(:,2),'.','color',colors{iExp})
    end
end
legend(expNames)

% discrimPairedBC = [pairedBC.e0; pairedBC.e5(:,:,1); pairedBC.e5(:,:,2)]; 
% allPairedBCN = [discrimPairedBC/max(discrimPairedBC(:)); pairedBC.e3/max(pairedBC.e3(:))];
% 
% [r p] = corr(allPairedBCN);

allPaired = [pairedBC.e0; pairedBC.e3; pairedBC.e5(:,:,1); pairedBC.e5(:,:,2)];
[r p] = corr(allPaired)

aa = [pairedBC.e0; pairedBC.e5(:,:,1); pairedBC.e5(:,:,2)];


figure
hold on
for iExp = [1 3]
    expName = expNames{iExp};
    for iNT = 1:size(pairedBC.(expName),3)
%         vals = pairedBC.(expName)(:,:,iNT)/max(discrimPairedBC(:));
        vals = pairedBC.(expName)(:,:,iNT);
        plot(vals(:,1),vals(:,2),'.','color',colors{iExp})
    end
end
% vals = pairedBC.e3/max(pairedBC.e3(:));
vals = pairedBC.e3;
plot(vals(:,1),vals(:,2),'.','color',colors{2})
xlabel('benefit')
ylabel('cost')

% colors are experiments, benefit targets are shapes, cost targets are shading
figure
hold on
for iExp = [1 3]
    expName = expNames{iExp};
    for iT = 1:size(pairedBC.(expName),1)
        for iNT = 1:size(pairedBC.(expName),3)
%             vals = pairedBC.(expName)(:,:,iNT)/max(discrimPairedBC(:));
            vals = pairedBC.(expName)(:,:,iNT);
            plot(vals(iT,1),vals(iT,2),shapes{iT},'color',colors{iExp},'MarkerFaceColor',faceColors{iNT})
        end
    end
end
% vals = pairedBC.e3/max(pairedBC.e3(:));
vals = pairedBC.e3;
for iT = 1:size(pairedBC.e3,1)
    plot(vals(iT,1),vals(iT,2),shapes{iT},'color',colors{2})
end
xlabel('benefit index')
ylabel('cost index')
xlim([-.4 1.2])
ylim([-.4 1.2])
set(gca,'XTick',[-.4 0 .4 .8 1.2])
set(gca,'YTick',[-.4 0 .4 .8 1.2])
axis square

% shapes are experiment, inner colors are benefit targets, outer colors are cost targets
figure
hold on
for iExp = [1 3]
    expName = expNames{iExp};
    tn = targetNames.(expName);
    alltargets = 1:numel(tn);
    for iT = 1:size(pairedBC.(expName),1)
        for iNT = 1:size(pairedBC.(expName),3)
            nontargets = setdiff(alltargets,iT);
            nt = nontargets(iNT);
%             vals = pairedBC.(expName)(:,:,iNT)/max(discrimPairedBC(:));
            vals = pairedBC.(expName)(:,:,iNT);
            plot(vals(iT,1),vals(iT,2),shapes{iExp},'color',colors{nt},'MarkerFaceColor',faceColors2{iT},'LineWidth',1.5,'MarkerSize',10)
        end
    end
end
% vals = pairedBC.e3/max(pairedBC.e3(:));
vals = pairedBC.e3;
for iT = 1:size(pairedBC.e3,1)
    nontargets = setdiff([1 2],iT);
    nt = nontargets;
    plot(vals(iT,1),vals(iT,2),shapes{2},'color',colors{nt},'MarkerFaceColor',faceColors2{iT},'LineWidth',1.5,'MarkerSize',10)
end
xlabel('benefit index')
ylabel('cost index')
xlim([-.4 1.2])
ylim([-.4 1.2])
set(gca,'XTick',[-.4 0 .4 .8 1.2])
set(gca,'YTick',[-.4 0 .4 .8 1.2])
axis square


% calculate correlation line and CIs
figure
x = allPaired(:,1); y = allPaired(:,2);
stats = regstats(y,x,'linear','beta');
% same as predint simultaneous functional confidence bounds in later matlab versions
[top, bot] = regression_line_ci(.05,stats.beta,x,y,100,min(x),max(x));

lims = [-.4 1.6];
% ticks = lims(1):.4:lims(end);
ticks = [0 .5 1 1.5];
figure 
hold on
plot(x,y,'.')
xx = linspace(min(x),max(x),101);
l = stats.beta(2)*xx + stats.beta(1);
plot(xx,l)
plot(xx,top)
plot(xx,bot)
plot(lims,lims,'k--')
xlim(lims)
ylim(lims)
set(gca,'XTick',ticks)
set(gca,'YTick',ticks)
axis square

% colors are experiments, targets not differentiated
% includes correlation line and CIs
figure
hold on
for iExp = 1:nExp
    expName = expNames{iExp};
    for iT = 1:size(pairedBC.(expName),1)
        for iNT = 1:size(pairedBC.(expName),3)
            vals = pairedBC.(expName)(:,:,iNT);
            plot(vals(iT,1),vals(iT,2),'.','color',colors{iExp})
        end
    end
end
plot(xx,l)
plot(xx,top)
plot(xx,bot)
plot(lims,lims,'k--')
xlabel('benefit index')
ylabel('cost index')
xlim(lims)
ylim(lims)
set(gca,'XTick',ticks)
set(gca,'YTick',ticks)
axis square

%% plot benefit/cost bar graphs for each experiment
for iExp = 1:nExp
    expName = expNames{iExp};
    tn = targetNames.(expName);
    nTargets = numel(tn);
    alltargets = 1:nTargets;

    pp = [];
    for iT = 1:nTargets
        vals = pairedBC.(expName);
        vals(:,2,:) = -1*vals(:,2,:); % make costs negative
        
        % order benefits and costs according to target number
        nontargets = setdiff(alltargets,iT);
        p = zeros(1,nTargets);
        p(iT) = vals(iT,1,1);
        p(nontargets) = vals(iT,2,:);
        pp(iT,:) = p;
    end

    figure
    bar(pp)
end
    
%% average costs across targets
if bootstrapErrorBars
    ciType = 68;
    switch ciType
        case 95
            ci = load('data/E0_resample_workspace_bootstrapPairedBC_20161128.mat','pairedBCCI');
            pairedBCCI.e0 = ci.pairedBCCI.e0;
            
            ci = load('data/E3_resample_workspace_bootstrapPairedBC_20161128.mat','pairedBCCI');
            pairedBCCI.e3 = ci.pairedBCCI.e3;
            
            ci = load('data/E5_resample_workspace_bootstrapPairedBC_20161128.mat','pairedBCAveCI');
            pairedBCCI.e5 = ci.pairedBCAveCI.e5;
        case 68
            ci = load('data/E0E3E5_resample_workspace_bootstrapPairedBC_CI68_20161128.mat','pairedBCCI','pairedBCAveCI');
            pairedBCCI.e0 = ci.pairedBCCI.e0;
            pairedBCCI.e3 = ci.pairedBCCI.e3;
            pairedBCCI.e5 = ci.pairedBCAveCI.e5;
    end
end

for iExp = 1:nExp
    expName = expNames{iExp};
    tn = targetNames.(expName);
    nTargets = numel(tn);
    alltargets = 1:nTargets;
    
    vals = pairedBC.(expName);
    vals(:,2,:) = -1*vals(:,2,:); % make costs negative
    
    bcave = mean(vals,3); % average costs across targets
    
    figure
    hold on
    h = bar(bcave);
    ylim([-1.2 1.2])
        
    if bootstrapErrorBars
        eb = pairedBCCI.(expName);
        eb(:,2,:) = -1*eb(:,2,:); % make costs negative
        
        for iT = 1:nTargets
            for iBC = 1:2
                xdata=get(get(h(iBC),'Children'),'XData');
                x=xdata(1,iT)+(xdata(3,iT)-xdata(1,iT))/2;
                y = bcave(iT,iBC);
                l = eb(iT,iBC,1)-y;
                u = eb(iT,iBC,2)-y;
                if iBC==2 % flip error bars for costs
                    ltemp = l;
                    l = u;
                    u = ltemp;
                end
                
                errorbar(x,y,l,u);
            end
        end
    end
end

%% average across targets plus costs for each target
for iExp = 1:nExp
    expName = expNames{iExp};
    tn = targetNames.(expName);
    nTargets = numel(tn);
    alltargets = 1:nTargets;
    
    vals = pairedBC.(expName);
    vals(:,2,:) = -1*vals(:,2,:); % make costs negative
    bcave = mean(vals,3); % average benefits and costs across targets
    
    % order benefits and costs according to target number
    if size(vals,3)==1
        pp = vals(:,:,1);
    else
        pp = [bcave vals(:,2,1) vals(:,2,2)];
    end

    figure
    bar(pp)
    ylim([-1 1])
end

%% stacked bars
for iExp = 1:nExp
    expName = expNames{iExp};
    vals = stackedBC.(expName);
    if strcmp(expName,'e5')
        sets = {[1 2],[1 3],[2 3]};
    else
        sets = {[1 2]};
    end
    
    figure
    for ip = 1:numel(sets)
        s = sets{ip};
        v(1,:) = vals(s(1),:,s(2));
        v(2,:) = vals(s(2),:,s(1));
        subplot(1,numel(sets),ip)
        bar(fliplr(v),'stacked')
        ylim([-.3 1.3])
    end
end
end

%% Randomization stats
if doRandomizationStats
    R = load('data/E0E3E5_resample_workspace_permutationPairedBC_20161128.mat');
    
    for iExp = 1:nExp
        expName = expNames{iExp};
        tn = targetNames.(expName);
        nTargets = numel(tn);
        for iT = 1:nTargets
            for iBC = 1:2
                if strcmp(expName,'e5')
                    val = mean(pairedBC.(expName)(iT,iBC,:),3); % mean across non-targets
                    nullDist = squeeze(R.pairedBCAveSamples.(expName)(iT,iBC,:));
                else
                    val = pairedBC.(expName)(iT,iBC);
                    nullDist = squeeze(R.pairedBCSamples.(expName)(iT,iBC,:));
                end
                pval.(expName)(iT,iBC) = nnz(abs(nullDist)>abs(val))/numel(nullDist);
            end
        end
    end
end

