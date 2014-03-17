function p = temporalAttentionParams

p.testingLocation = 'CarrascoL1'; % 'CarrascoL1','laptop'

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
        p.screenRes = [1600 1200];
        p.viewDist = 57;
    otherwise
        error('Testing location not found in temporalAttentionParams.')
end     

p.keyCodes = KbName(p.keyNames);
p.backgroundColor = 0.5;
p.nReps = 3;
p.font = 'Verdana';

% Conditions
p.targetContrasts = [.16 .32];
p.respInterval = [1 2]; % [1=early 2=late]
p.cueValidity = [1 -1 0]; % [1=valid -1=invalid 0=neutral]
% p.propValid = 0.67;
% p.cueValidityFactor = generatePropFactor(p.propValid);
p.cueValidityFactor = [1 1 1 2]; % 50% valid, 25% invalid, 25% neutral
p.propValid = nnz(p.cueValidityFactor==1)./nnz(p.cueValidityFactor<3);
p.propNeutral = nnz(p.cueValidityFactor==3)./numel(p.cueValidityFactor);

% Timing
p.soas = [300 733]/1000; % [short long]
p.preCueDur = 0.5; % time between fixation onset and cue
p.cueDur = 0.2;
p.targetDur = 2/60; % 33 ms
p.maskSOA = 4/60; % time between target onset and mask onset
p.maskDur = 1/60;
p.respCueSOA = p.soas(2) + 0.5;
p.iti = 0.5; % inter-trial interval

% Images
p.imPos = [4 4];
p.imSize = [4 4]; % this is the size of the image container that holds the stim
p.targetSize = 0.5; % sigma of gaussian
p.spatialFrequency = 4;
p.targetOrientations = [-10 10];

p.maskType = 'filterednoise'; % whitenoise, verticalgrating, crossedgratings, filterednoise
p.maskContrast = 1;

% Sounds
p.Fs = 44100;
p.cueFreqs = [784 523]; % [higher G = target 1, lower C = target 2]
for iTone = 1:numel(p.cueFreqs)
    tone = MakeBeep(p.cueFreqs(iTone), p.cueDur, p.Fs);
    p.cueTones(iTone,:) = applyEnvelope(tone, p.Fs);
end
% 10^0.5 for every 10dB



