function praise = highpraise(idx)
%
% function praise = highpraise([idx])

% Note, weights are a little weird. This might be true sometimes:
% .5 + .5*(weight/sum(weights)) = trueprop;

if nargin==0
    idx = [];
end

% praises = {1, 2, 3};
praises = {'Great work!','Rock star!','Champ!'};
weights = [4 1 1];

if isempty(idx)
    x = rand(1,numel(praises)).*weights;
    choice = x==max(x);
else
    choice = idx;
end

praise = praises{choice};


% % test weights
% for i=1:1000
% p(i) = highpraise;
% end