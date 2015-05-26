function [f xgrid] = DoG(mu, sd, amp, b, xgrid)

if nargin < 1 || isempty(mu)
    mu = 0;
end
if nargin < 2 || isempty(sd)
    sd = 1;
end
if nargin < 3 || isempty(amp)
    amp = 1;
end
if nargin < 4 || isempty(b)
    b = 0;
end
if nargin < 5 || isempty(xgrid)
    xgrid = -10:.01:10;
end

g = normpdf(xgrid, mu, sd);
dg = diff(g);
f = dg./max(dg).*amp + b;

f = [f 0]; % restore to original length