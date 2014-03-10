function p = temporalAttentionParams

p.testingLocation = 'laptop'; % 'CarrascoL1','laptop'

switch p.testingLocation
    case 'laptop'
        p.keyNames = {'1!','2@'};
        p.refRate = 1/60;
        p.screenSize = [9 13]; % (in)
        p.screenRes = [900 1440];
        p.viewDist = 21; % (in)
    case 'CarrascoL1'
        p.keyNames = {'1!','2@'};
        p.refRate = 1/60;
        p.screenSize = [40 30];
        p.screenRes = [1024 768];
        p.viewDist = 57;
    otherwise
        error('Testing location not found in temporalAttentionParams.')
end     

p.keyCodes = KbName(p.keyNames);
p.backgroundColor = 0.5;
p.nReps = 1;

% Conditions
p.targetContrasts = [0.16 0.32];
p.cuedInterval = [1 2]; % [1=early 2=late]
p.cueValidity = [1 0]; % [1=valid 0=invalid]
p.propValid = 0.75;
p.cueValidityFactor = generatePropFactor(p.propValid);

% Timing
p.soas = [300 550]/1000; % [short long]
p.targetDur = 2/60; % 33 ms
p.respCueSOA = p.soas(2) + 0.3;

% Images
p.fixSize = 0.1;
p.imPos = [4 4];
p.imSize = [6 6]; % this is the size of the image container that holds the stim
p.targetSize = 0.5; % sigma of gaussian
p.spatialFrequency = 4;
p.targetOrientations = [-15 15];

% Sounds
v = 1:1000;
envelope = [0:.05:1 ones(1,1000-42) 1:-.05:0];
p.startSound = (0.25 * sin(3*pi*v/30)).*envelope;






