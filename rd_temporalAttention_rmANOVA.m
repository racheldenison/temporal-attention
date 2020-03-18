D = load('data/E2_SOA_cbD6_run98_N5_workspace_20191003.mat');

% data = D.dpData;
data = D.rtData;

allData(:,:,:,1) = shiftdim(data{1},2);
allData(:,:,:,2) = shiftdim(data{2},2);

[tbl, rm] = simple_mixed_anova(allData,[],{'validity','soa','target'});

% t = mauchly(rm)

t1bData = allData(:,[1 3],:,1);
t1cData = allData(:,[2 3],:,1);
t2bData = allData(:,[1 3],:,2);
t2cData = allData(:,[2 3],:,2);

[tbl, rm] = simple_mixed_anova(t1bData,[],{'validity','soa'});
[tbl, rm] = simple_mixed_anova(t2cData,[],{'validity','soa'});


viDiffData = squeeze(allData(:,1,:,:) - allData(:,2,:,:));

[tbl, rm] = simple_mixed_anova(viDiffData,[],{'soa','target'});