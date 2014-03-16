function rd_saveAllFigs(f, figNames, figNamePrefix, figDir, fileType)
%
% function rd_saveAllFigs(h, figNames, [figNamePrefix], [figDir], [fileType])
%
% h is a vector of figure handles
% figNames is a cell array of figure names
% figNamePrefix is a string that will be the start of all figure names
% (default is no prefix)
% figDir is the figure directory path (default is 'figures')
% fileType is a string giving the 'device' (default is '-dpng')
%
% Rachel Denison
% 2014 March 14 (pi day!)

if nargin < 3 || isempty(figNamePrefix)
    figNamePrefix = '';
end
if nargin < 4 || isempty(figDir)
    figDir = 'figures';
end
if nargin < 5 || isempty(fileType)
    fileType = '-dpng';
end

for iF = 1:numel(f)
    figFile = sprintf('%s/%s_%s', figDir, figNamePrefix, figNames{iF});
    print(f(iF), fileType, figFile)
end