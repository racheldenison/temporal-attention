function [dprime, criterion] = getDprimeFromTrials(trials, headers)

targetTypeIdx = strcmp(headers, 'respTargetState');
accIdx = strcmp(headers, 'correct');

targetType = sign(trials(:,targetTypeIdx));
acc = trials(:,accIdx);

targetTypes = unique(targetType);

if numel(targetTypes)~=2
    error('There should be 2 target types.')
end

% let target type 1 be the "signal" and type 2 be the "noise"
signal = targetType==targetTypes(1);
noise = targetType==targetTypes(2);

if nnz(signal)~=nnz(noise)
    error('We expect equal numbers of signal and noise trials.')
end

hit = signal & acc==1;
fa = noise & acc==0;

[dprime, criterion] = rd_dprime2(nnz(hit),nnz(fa),nnz(signal),nnz(noise));