% rd_plotTemporalAttentionTradeoffsMean.m

%% setup
e0 = load('data/E0_workspace_run09_N10_20160224.mat');
e3 = load('data/E3_workspace_run09_N12_20160224.mat');
e5 = load('data/E5_workspace_run01_N10_20160806.mat');

% fractional benefits and costs?
normalizeWithinTarget = 1;
baselineOpt = 'vi'; % 'vi','condsum'

% valid vs. invalid
vi.e0 = squeeze(mean(e0.accDataCP(1,:,:),2));
vi.e3 = squeeze(mean(e3.pdData(1,:,:),2));
vi.e5 = squeeze(mean(e5.accDataCP(1,:,:),2));

vi12.e5 = squeeze(mean(e5.accDataIBP(1:2,:,:),2))'; % target x VI1/VI2

% sum across conditions
condsum.e0 = sum(e0.accMeanC)';
condsum.e3 = sum(e3.paramsMean.sd)';
condsum.e5 = sum(e5.accMeanC)';

condsum12.e5(:,1) = sum(e5.accMeanIB(1:3,:)); % I1
condsum12.e5(:,2) = sum(e5.accMeanIB([1 2 4],:)); % I2

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
        bn = sprintf('%sb', tn{iT});
        benefit = bc.(expName).(bn);
        if normalizeWithinTarget
            benefit = benefit/base.(expName)(iT);
        end
        nontargets = setdiff(alltargets,iT);
        for iNT = 1:numel(nontargets)
            nt = nontargets(iNT);
            cn = sprintf('%sc_cue%s',tn{nt},tn{iT});
            cost = bc.(expName).(cn);
            if normalizeWithinTarget
                if strcmp(expName,'e5')
                    targets = setdiff(alltargets,nt);
                    idx = find(targets==iT);
                    cost = cost/base12.(expName)(nt,idx);
                else
                    cost = cost/base.(expName)(nt);
                end
            end
            pairedBC.(expName)(iT,:,iNT) = [benefit cost];
            pairNames.(expName){iT,iNT} = sprintf('%s_%s',bn,cn);
        end
    end
end

%% plot
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

% colors are experiments, targets not differentiated
figure
hold on
for iExp = [1 3]
    expName = expNames{iExp};
    for iT = 1:size(pairedBC.(expName),1)
        for iNT = 1:size(pairedBC.(expName),3)
            vals = pairedBC.(expName)(:,:,iNT);
            plot(vals(iT,1),vals(iT,2),'.','color',colors{iExp})
        end
    end
end
vals = pairedBC.e3;
for iT = 1:size(pairedBC.e3,1)
    plot(vals(iT,1),vals(iT,2),'.','color',colors{2})
end
xlabel('benefit index')
ylabel('cost index')
xlim([-.4 1.2])
ylim([-.4 1.2])
set(gca,'XTick',[-.4 0 .4 .8 1.2])
set(gca,'YTick',[-.4 0 .4 .8 1.2])
axis square

x = allPaired(:,1); y = allPaired(:,2);
stats = regstats(y,x,'linear','beta');
% same as predint simultaneous functional confidence bounds in later matlab versions
[top, bot] = regression_line_ci(.05,stats.beta,x,y,100,-.4,1.2);

figure 
hold on
plot(x,y,'.')
xx = linspace(-.4,1.2,101);
l = stats.beta(2)*xx + stats.beta(1);
plot(xx,l)
plot(xx,top)
plot(xx,bot)
plot(xx,xx,'k--')
xlim([-.4 1.2])
ylim([-.4 1.2])
set(gca,'XTick',[-.4 0 .4 .8 1.2])
set(gca,'YTick',[-.4 0 .4 .8 1.2])
axis square


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



