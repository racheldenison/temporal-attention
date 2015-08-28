function turnwhite(f)
%
% turnwhite([f])
%
% turns the background color of a figure white
% f is the figure handle
% if f is not given, use the current figure

if nargin==0
    f = gcf;
end

set(f,'color','w')