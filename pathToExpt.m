function exptpath = pathToExpt(directory)

% exptpath = pathToExpt(directory)

exptpath = '~/Documents/NYU/Grants_&_Apps/NSF_SBE_Postdoc_Fellowship';
% exptpath = sprintf('%s/Temporal_Attention', pathToCarrascoExpts);

if nargin==1
    exptpath = sprintf('%s/%s', exptpath, directory);
end