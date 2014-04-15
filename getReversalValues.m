function rvals = getReversalValues(vals)

% function rvals = getReversalValues(vals)
%
% Given a series of values, find 
%
% Rachel Denison
% 2014 April 15

% if size(vals,2)>size(vals,1)
%     vals = vals'; % doesn't matter if 1D
% end
% 1st diff
a = diff(vals);
vals = vals(2:end);
% get rid of zeros
b = a;
b(a==0) = [];
vals(a==0) = [];
% 2nd diff
c = diff(b);
vals = vals(2:end);
% get reversal values
rvals = vals(c~=0);