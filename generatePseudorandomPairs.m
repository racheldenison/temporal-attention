function prps = generatePseudorandomPairs(nLevels, nPairs, noPairEquality)

% function prps = generatePseudorandomPairs(nLevels, nPairs,
% [noPairEquality=0])
%
% Creates a randomized [nPairs x 2] list of pairs with values taking 
% 1:nLevels. Attempts roughly equal numbers of each pair. 
% If noPairEquality = 1, then the two members of a pair must be
% different.
%
% Rachel Denison
% 14 April 2014

if nargin < 3
    noPairEquality = 0;
end

pairs = fullfact([nLevels nLevels]);
if noPairEquality
    same = pairs(:,1)==pairs(:,2);
    pairs(same,:) = [];
end
prps0 = repmat(pairs,ceil(nPairs/size(pairs,1)),1);
randOrder = randperm(size(prps0,1));
prps = prps0(randOrder(1:nPairs),:);
