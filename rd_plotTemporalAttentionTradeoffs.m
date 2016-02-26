% rd_plotTemporalAttentionTradeoffs.m

%% setup
e0 = load('data/E0_workspace_run09_N10_20160224.mat');
e3 = load('data/E3_workspace_run09_N12_20160224.mat');

pdata.e0 = e0.accDataCP;
pdata.e3 = e3.pdData;

expNames = fields(pdata);
nExp = numel(expNames);

%% normalize
for iE = 1:nExp
    exp = expNames{iE};
    pdatan.(exp) = zscore(pdata.(exp),0,2);
end

pdatanAll = cat(2, pdatan.e0, pdatan.e3);

[rT1bT2c pT1bT2c] = corr(pdatanAll(2,:,1)', pdatanAll(3,:,2)');
[rT1cT2b pT1cT2b] = corr(pdatanAll(3,:,1)', pdatanAll(2,:,2)');

%% benefit-cost index
for iE = 1:nExp
    exp = expNames{iE};
    bcIndex.(exp)(1,:) = pdata.(exp)(2,:,1) - pdata.(exp)(3,:,1); % T1 B-C
    bcIndex.(exp)(2,:) = pdata.(exp)(3,:,2) - pdata.(exp)(2,:,2); % T2 C-B
end

% normalize
for iE = 1:nExp
    exp = expNames{iE};
    bcIndexN.(exp) = zscore(bcIndex.(exp),0,2);
end

bcIndexNAll = cat(2, bcIndexN.e0, bcIndexN.e3);

[rBCIndex pBCIndex] = corr(bcIndexNAll(1,:)', bcIndexNAll(2,:)');
    
%% plot
colors = {'b','r'};

%% T1 benefit vs. T2 cost
figure
hold on
for iE = 1:nExp
    subplot(1,nExp,iE)
    hold on
    exp = expNames{iE};
    plot(pdata.(exp)(2,:,1), pdata.(exp)(3,:,2),'.','color', colors{iE}); % T1 benefit vs. T2 cost
    xlabel('T1 benefit')
    ylabel('T2 cost')
    title(exp)
    switch exp
        case 'e0'
            plot([-.2 .2],[0 0],'--k')
            plot([0 0],[-.2 .2],'--k')
            xlim([-.2 .2])
            ylim([-.2 .2])
        case 'e3'
            plot([-15 15],[0 0],'--k')
            plot([0 0],[-15 15],'--k')
            xlim([-15 15])
            ylim([-15 15])
    end
end

figure
hold on
for iE = 1:nExp
    exp = expNames{iE};
    plot(pdatan.(exp)(2,:,1), pdatan.(exp)(3,:,2),'.','color', colors{iE}); % T1 benefit vs. T2 cost
end
xlim([-3 3])
ylim([-3 3])
legend(expNames)
title('T1 benefit vs. T2 cost, normalized')

figure
plot(pdatanAll(2,:,1), pdatanAll(3,:,2),'.'); % T1 benefit vs. T2 cost
title('T1 benefit vs. T2 cost, normalized')

%% T1 cost vs. T2 benefit
figure
hold on
for iE = 1:nExp
    subplot(1,nExp,iE)
    hold on
    exp = expNames{iE};
    plot(pdata.(exp)(3,:,1), pdata.(exp)(2,:,2),'.','color', colors{iE}); % T1 benefit vs. T2 cost
    xlabel('T1 cost')
    ylabel('T2 benefit')
    title(exp)
    switch exp
        case 'e0'
            plot([-.2 .2],[0 0],'--k')
            plot([0 0],[-.2 .2],'--k')
            xlim([-.2 .2])
            ylim([-.2 .2])
        case 'e3'
            plot([-15 15],[0 0],'--k')
            plot([0 0],[-15 15],'--k')
            xlim([-15 15])
            ylim([-15 15])
    end
end

figure
hold on
for iE = 1:nExp
    exp = expNames{iE};
    plot(pdatan.(exp)(3,:,1), pdatan.(exp)(2,:,2),'.','color', colors{iE}); % T1 benefit vs. T2 cost
end
xlim([-3 3])
ylim([-3 3])
legend(expNames)
title('T1 cost vs. T2 benefit, normalized')

figure
plot(pdatanAll(3,:,1), pdatanAll(2,:,2),'.'); % T1 cost vs. T2 benefit
title('T1 cost vs. T2 benefit, normalized')

%% BC Index
figure
hold on
for iE = 1:nExp
    subplot(1,nExp,iE)
    exp = expNames{iE};
    plot(bcIndex.(exp)(1,:), bcIndex.(exp)(2,:),'.','color', colors{iE}); % T1 index vs. T2 index
    xlabel('T1 BC index')
    ylabel('T2 BC index')
    title(exp)
end

figure
hold on
for iE = 1:nExp
    exp = expNames{iE};
    plot(bcIndexN.(exp)(1,:), bcIndexN.(exp)(2,:),'.','color', colors{iE});
end
legend(expNames)
title('T1 vs. T2 BC index, normalized')

figure
plot(bcIndexNAll(1,:), bcIndexNAll(2,:),'.'); 
title('T1 vs. T2 BC index, normalized')


