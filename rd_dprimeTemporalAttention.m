function [dprime, criterion] = rd_dprimeTemporalAttention(trials, trials_headers)

%% setup
stateIdx = find(strcmp(trials_headers,'respTargetState'));
validityIdx = find(strcmp(trials_headers,'cueValidity'));
targetIdx = find(strcmp(trials_headers,'respInterval'));
accIdx = find(strcmp(trials_headers,'correct'));

targetState = sign(trials(:,stateIdx));
states = unique(targetState);

cueValidity = trials(:,validityIdx);
cvs = unique(cueValidity);

target = trials(:,targetIdx);
targets = unique(target);

acc = trials(:,accIdx);

%% proportion correct for each target state & condition
for iT = 1:numel(targets)
    for iV = 1:numel(cvs)
        for iS = 1:numel(states)
            w = target == targets(iT) & ...
                cueValidity == cvs(iV) & ...
                targetState == states(iS);
            pc(iV,iS,iT) = mean(acc(w));
            nt(iV,iS,iT) = nnz(w);
        end
    end
end

%% dprime and criterion
% equivalent to norminv(pc(:,1)) + norminv(pc(:,2))
h = squeeze(pc(:,1,:));
fa = squeeze(1-pc(:,2,:));
nt = squeeze(nt(:,1,:));

% need h < 1, f > 0
% change the proportion by the magnitude of one trial if needed
h(h==1) = 1 - 1./(nt(h==1));
fa(fa==0) = 1./(nt(fa==0));

zh = norminv(h); 
zfa = norminv(fa);

dprime = zh - zfa;
criterion = -0.5*(zh+zfa);

