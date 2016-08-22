% rd_plotTemporalAttentionTradeoffsMean.m

%% setup
e0 = load('data/E0_workspace_run09_N10_20160224.mat');
e3 = load('data/E3_workspace_run09_N12_20160224.mat');
e5 = load('data/E5_workspace_run01_N10_20160806.mat');

% benefits and costs
% VI, VN, NI
bc.e0.t1b = e0.accDataCP(2,:,1);
bc.e0.t1c_cuet2 = e0.accDataCP(3,:,1);
bc.e0.t2b = e0.accDataCP(2,:,2);
bc.e0.t2c_cuet1 = e0.accDataCP(3,:,2);

bc.e3.t1b = e3.pdData(2,:,1);
bc.e3.t1c_cuet2 = e3.pdData(3,:,1);
bc.e3.t2b = e3.pdData(2,:,2);
bc.e3.t2c_cuet1 = e3.pdData(3,:,2);

% VI1, VI2, VN, NI1, NI2
bc.e5.t1b = e5.accDataIBP(3,:,1);
bc.e5.t2b = e5.accDataIBP(3,:,2);
bc.e5.t3b = e5.accDataIBP(3,:,3);
bc.e5.t1c_cuet2 = e5.accDataIBP(4,:,1);
bc.e5.t1c_cuet3 = e5.accDataIBP(5,:,1);
bc.e5.t2c_cuet1 = e5.accDataIBP(4,:,2);
bc.e5.t2c_cuet3 = e5.accDataIBP(5,:,2);
bc.e5.t3c_cuet1 = e5.accDataIBP(4,:,3);
bc.e5.t3c_cuet2 = e5.accDataIBP(5,:,3);

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
        nontargets = setdiff(alltargets,iT);
        for iNT = 1:numel(nontargets)
            nt = nontargets(iNT);
            cn = sprintf('%sc_cue%s',tn{nt},tn{iT});
            cost = bc.(expName).(cn);
            pairedBC.(expName)(:,:,iT,iNT) = [benefit' cost'];
            pairNames.(expName){iT,iNT} = sprintf('%s_%s',bn,cn);
        end
    end
end

%% plot
colors = {'b','r','k'};
% colors = {[97 47 255]/255,[255 4 0]/255,[255 177 6]/255};
% colors = {[.1 .1 .1], [.4 .4 .4], [.7 .7 .7]};
shapes = {'o','s','^'};
faceColors = {[1 1 1],[.5 .5 .5]};
faceColors2 = colors;
figure
hold on
allVals.all = [];
for iExp = 1:nExp
%     figure
%     hold on
    expName = expNames{iExp};
    allVals.(expName) = [];
    dataMax = max(pairedBC.(expName)(:));
%     dataMean = mean(pairedBC.(expName)(:));
%     dataSD = std(pairedBC.(expName)(:));
    for iT = 1:size(pairedBC.(expName),3)
%     for iT = 2
        for iNT = 1:size(pairedBC.(expName),4)
%         for iNT = 1
            vals = squeeze(pairedBC.(expName)(:,:,iT,iNT));
%                     vals = pairedBC.(expName)(:,:,iNT)/dataMax;
            %         vals = (pairedBC.(expName)(:,:,iNT)-dataMean)/dataSD;
            plot(vals(:,1),vals(:,2),'.','color',colors{iExp})
            allVals.all = [allVals.all; vals];
            allVals.(expName) = [allVals.(expName); vals];
        end
    end
end
legend(expNames)

discrimPairedBC = [allVals.e0; allVals.e5]; 
allPairedBCN = [discrimPairedBC/max(discrimPairedBC(:)); allVals.e3/max(allVals.e3(:))];

[r p] = corr(allPairedBCN);

figure
hold on
for iExp = [1 3]
    expName = expNames{iExp};
    for iT = 1:size(pairedBC.(expName),3)
        for iNT = 1:size(pairedBC.(expName),4)
            vals = pairedBC.(expName)(:,:,iT,iNT)/max(discrimPairedBC(:));
            plot(vals(:,1),vals(:,2),'.','color',colors{iExp})
        end
    end
end
for iT = 1:size(pairedBC.e3,3)
    vals = pairedBC.e3(:,:,iT)/max(pairedBC.e3(:));
    plot(vals(:,1),vals(:,2),'.','color',colors{2})
end
xlabel('benefit')
ylabel('cost')

figure
hold on
for iExp = [1 3]
    expName = expNames{iExp};
    for iT = 1:size(pairedBC.(expName),3)
        for iNT = 1:size(pairedBC.(expName),4)
            vals = pairedBC.(expName)(:,:,iT,iNT)/max(discrimPairedBC(:));
            plot(vals(:,1),vals(:,2),shapes{iT},'color',colors{iExp},'MarkerFaceColor',faceColors{iNT})
        end
    end
end
for iT = 1:size(pairedBC.e3,3)
    vals = pairedBC.e3(:,:,iT)/max(pairedBC.e3(:));
    plot(vals(:,1),vals(:,2),shapes{iT},'color',colors{2})
end
xlabel('benefit')
ylabel('cost')
xlim([-.8 1])
ylim([-.8 1])
axis square

% figure
% hold on
% for iExp = [1 3]
%     expName = expNames{iExp};
%     tn = targetNames.(expName);
%     alltargets = 1:numel(tn);
%     for iT = 1:size(pairedBC.(expName),1)
%         for iNT = 1:size(pairedBC.(expName),3)
%             nontargets = setdiff(alltargets,iT);
%             nt = nontargets(iNT);
%             vals = pairedBC.(expName)(:,:,iNT)/max(discrimPairedBC(:));
%             plot(vals(iT,1),vals(iT,2),shapes{iExp},'color',colors{iT},'MarkerFaceColor',faceColors2{nt},'LineWidth',1.5,'MarkerSize',10)
%         end
%     end
% end
% vals = pairedBC.e3/max(pairedBC.e3(:));
% for iT = 1:size(pairedBC.e3,1)
%     nontargets = setdiff([1 2],iT);
%     nt = nontargets;
%     plot(vals(iT,1),vals(iT,2),shapes{2},'color',colors{iT},'MarkerFaceColor',faceColors2{nt},'LineWidth',1.5,'MarkerSize',10)
% end
% xlabel('benefit')
% ylabel('cost')
% xlim([-.2 1.01])
% ylim([-.2 1.01])
% axis square



