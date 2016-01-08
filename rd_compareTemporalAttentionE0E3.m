% rd_compareTemporalAttentionE0E3.m

%% setup
e0 = load('data/E0_workspace_run09_N10_20151210.mat');
e3 = load('data/adjust_fit_group_stats_mixtureNoBiasMaxPosterior_run09_N12_20150512.mat');

e3.subjectInits = {'bl','rd','id','ec','ld','en','sj','ml','ca','jl','ew','jx'};

subjects = intersect(e0.subjectInits, e3.subjectInits);
nSubjects = numel(subjects);

%% measures
measureNameE0 = 'accDataC';
measureNameE3 = 'g';

measureE0 = e0.(measureNameE0);
measureE3 = e3.paramsData.(measureNameE3);

%% get values for common subjects
for iS = 1:nSubjects
    subject = subjects{iS};
    for iT = 1:2
        valsE0(:,iT,iS) = measureE0{iT}(:,strcmp(e0.subjectInits, subject));
        valsE3(:,iT,iS) = measureE3(:,iT,strcmp(e3.subjectInits, subject));
    end
end

%% calculate differences between conditions
% valid vs. invalid
viE0 = squeeze(valsE0(1,:,:) - valsE0(2,:,:))';
viE3 = squeeze(valsE3(2,:,:) - valsE3(1,:,:))';

% valid vs. neutral
vnE0 = squeeze(valsE0(1,:,:) - valsE0(3,:,:))';
vnE3 = squeeze(valsE3(3,:,:) - valsE3(1,:,:))';

% valid vs. invalid
niE0 = squeeze(valsE0(3,:,:) - valsE0(2,:,:))';
niE3 = squeeze(valsE3(2,:,:) - valsE3(3,:,:))';

%% plot figs
figure
subplot(1,3,1)
plot(viE0, viE3, '.')
xlabel(measureNameE0)
ylabel(measureNameE3)
title('valid vs. invalid')
legend('T1','T2','location','best')

subplot(1,3,2)
plot(vnE0, vnE3, '.')
title('valid vs. neutral')

subplot(1,3,3)
plot(niE0, niE3, '.')
title('neutral vs. invalid')


