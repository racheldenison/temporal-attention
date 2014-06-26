function exptpath = pathToExpt(directory)

% exptpath = pathToExpt(directory)

exptpath = sprintf('%s/Temporal_Attention', pathToCarrascoExpts);

if nargin==1
    exptpath = sprintf('%s/%s', exptpath, directory);
end