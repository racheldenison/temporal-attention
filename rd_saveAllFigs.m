function rd_saveAllFigs(f, figNames, figNamePrefix, figDir, fileType)
%
% function rd_saveAllFigs(f, figNames, [figNamePrefix], [figDir], [fileType])
%
% f is a vector of figure handles. [] to save all open figures.
% figNames is a cell array of figure names
% figNamePrefix is a string that will be the start of all figure names
% (default is no prefix)
% figDir is the figure directory path (default is 'figures')
% fileType is a string giving the 'device' (default is '-dpng')
%
% Rachel Denison
% 2014 March 14 (pi day!)

if nargin < 3 || isempty(figNamePrefix)
    prefix = '';
else
    prefix = sprintf('%s_', figNamePrefix);
end
if nargin < 4 || isempty(figDir)
    figDir = 'figures';
end
if nargin < 5 || isempty(fileType)
    fileType = '-dpng';
end

if isempty(f)
    f = sort(findobj('Type','figure'));
end
if numel(f)~=numel(figNames)
    fprintf('\nSorry, the number of open figures does not match the number of figure names. Not saving anything.\n\n')
    return
end

for iF = 1:numel(f)
    figFile = sprintf('%s/%s%s', figDir, prefix, figNames{iF});
    print(f(iF), fileType, '-r0', figFile)
end