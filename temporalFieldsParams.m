function p = temporalFieldsParams

p.testingLocation = 'CarrascoL1'; % 'CarrascoL1','laptop'

switch p.testingLocation
    case 'laptop'
        p.keyNames = {'1!','2@'};
    case 'CarrascoL1'
        p.keyNames = {'1!','2@'};
        p.refRate = 1/120;
        p.screenSize = [40 30];
        p.screenRes = [1024 768];
        p.viewDist = 57;
    otherwise
        error('Testing location not found in temporalFieldsParams.')
end     

p.keyCodes = KbName(p.keyNames);
p.backgroundColor = 0.5;
p.nReps = 5;

% SOA conditions
p.targOnlyCode = 0.5; % SOA code for target only condition
p.soas = [-20 -10:2:-2 2:2:10 20]*1/60;
p.soas = [p.soas p.targOnlyCode]; % p.targOnlyCode means include surround absent condition
% p.soas = [p.targOnlyCode p.targOnlyCode]; 

% Target present / absent conditions
p.targetPresAbs = [1 0]; % 1=present, 0=absent

% Stimulus params
p.imageDur = 1/60;
p.cueTargetSOA = 0.5;
p.postStimCushion = 0.2;

p.imPos = [2 -2];
p.imSize = [4 4]; % this is the size of the image container that holds the stim
p.fixSize = 0.1;
p.targetSize = 2;
p.surroundSize = 3;
p.spatialFrequency = 3;
p.orientation = 0;
p.targetContrast = 0.05;
p.surroundContrast = 0.1;





